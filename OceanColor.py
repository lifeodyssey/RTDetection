# DataProducts/NASA/OceanColor.py - DatasetCollections and Grids that
# represent products published by the NASA GSFC OceanColor Group
#
# Copyright (C) 2010 Jason J. Roberts
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License (available in the file LICENSE.TXT)
# for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

import datetime
import os
import platform
import sys
import time 

from GeoEco.Datasets import Dataset, QueryableAttribute, DatasetCollection, Grid
from GeoEco.Datasets.Virtual import TimeSeriesGridStack, GridSliceCollection, RotatedGlobalGrid, ClippedGrid, ClimatologicalGridCollection
from GeoEco.Dependencies import PythonAggregatedModuleDependency, SoftwareNotInstalledError
from GeoEco.DynamicDocString import DynamicDocString
from GeoEco.Internationalization import _, UnicodeToUserPreferredEncoding, UserPreferredEncodingToUnicode
from GeoEco.Logging import Logger


class OceanColorLevel3SMIFileSearcher(DatasetCollection):
    __doc__ = DynamicDocString()

    def __init__(self, sensor=None, temporalResolution=None, product=None, startDate=None, endDate=None, timeout=60, maxRetryTime=300, cacheDirectory=None):
        self.__doc__.Obj.ValidateMethodInvocation()

        # Initialize our properties.

        self._Sensor = sensor
        self._TemporalResolution = temporalResolution
        self._Product = product
        self._StartDate = startDate
        self._EndDate = endDate
        self._Timeout = timeout
        self._MaxRetryTime = maxRetryTime

        self._NetCDF4Available = None
        self._QueryResults = None
        self._Http = None

        # Define a function used to calculate the value of the EndDate
        # queryable attribute from the value of the DateTime queryable
        # attribute.
        
        def _GetEndDate(selfOrDict, startDate):

            # Get the temporal resolution.
            
            temporalResolution = None
            if isinstance(selfOrDict, dict):
                if u'TemporalResolution' in selfOrDict:
                    temporalResolution = selfOrDict[u'TemporalResolution']
            else:
                temporalResolution = selfOrDict.GetQueryableAttributeValue(u'TemporalResolution')

            # Calculate the end date string from the temporal resolution.

            if temporalResolution is not None:
                temporalResolution = temporalResolution.lower()
                
                if temporalResolution == u'daily':
                    return ''
                
                if temporalResolution == u'8day':
                    return unicode((startDate + datetime.timedelta(days=7)).strftime('%Y%j'))
                
                if temporalResolution == u'monthly':
                    if startDate.month == 12:
                        return unicode((datetime.datetime(startDate.year + 1, 1, 1) - datetime.timedelta(days=1)).strftime('%Y%j'))
                    return unicode((datetime.datetime(startDate.year, startDate.month + 1, 1) - datetime.timedelta(days=1)).strftime('%Y%j'))
                
                if temporalResolution == u'annual':
                    return unicode((datetime.datetime(startDate.year + 1, 1, 1) - datetime.timedelta(days=1)).strftime('%Y%j'))
                    
            return None

        # Define the allowed queryable attributes.

        queryableAttributes = (QueryableAttribute(u'Sensor', _(u'Sensor'), UnicodeStringTypeMetadata(allowedValues=[u'aqua', u'czcs', u'meris', u'octs', u'seawifs', u'terra'], makeLowercase=True)),
                               QueryableAttribute(u'SensorCode', _(u'Sensor abbreviation code'), UnicodeStringTypeMetadata(allowedValues=[u'a', u'c', u'o', u's', u't'], makeLowercase=True)),
                               QueryableAttribute(u'TemporalResolution', _(u'Temporal resolution'), UnicodeStringTypeMetadata(allowedValues=[u'daily', u'8day', u'monthly', u'annual'], makeLowercase=True)),
                               QueryableAttribute(u'TemporalResolutionCode', _(u'Temporal resolution abbreviation code'), UnicodeStringTypeMetadata(allowedValues=[u'DAY', u'8D', u'MO', u'YR'], makeLowercase=True)),
                               QueryableAttribute(u'SpatialResolution', _(u'Spatial resolution'), UnicodeStringTypeMetadata(allowedValues=[u'4km', u'9km'], makeLowercase=True)),
                               QueryableAttribute(u'SpatialResolutionCode', _(u'Spatial resolution abbreviation code'), UnicodeStringTypeMetadata(allowedValues=[u'4', u'9'], makeLowercase=True)),
                               QueryableAttribute(u'ProductCode', _(u'Level 3 SMI product code'), UnicodeStringTypeMetadata()),
                               QueryableAttribute(u'DateTime', _(u'Start date'), DateTimeTypeMetadata()),
                               QueryableAttribute(u'EndDate', _(u'End date string'), UnicodeStringTypeMetadata(mustMatchRegEx=ur'\d\d\d\d\d\d\d\d'), derivedFromAttr=u'DateTime', derivedValueFunc=_GetEndDate))

        # Initialize the base class.

        super(OceanColorLevel3SMIFileSearcher, self).__init__(queryableAttributes=queryableAttributes, cacheDirectory=cacheDirectory)

    def _Close(self):
        if hasattr(self, '_Http') and self._Http is not None:
            self._LogDebug(_(u'%(class)s 0x%(id)08X: Closed.'), {u'class': self.__class__.__name__, u'id': id(self)})
            self._Http = None
        super(OceanColorLevel3SMIFileSearcher, self)._Close()

    def _GetDisplayName(self):
        return _(u'NASA OceanColor data archive')

    def _QueryDatasets(self, parsedExpression, progressReporter, options, parentAttrValues):

        # If we have not queried the server yet, do it now.

        if self._QueryResults is None:
            self._QueryServer()

        # Go through the list of URLs returned by the server, testing
        # whether each one matches the query expression. For each
        # match, construct a NetCDFFile or HDF4SDSCollection instance,
        # query it, and return the resulting Grids.

        from GeoEco.Datasets.NetCDF import NetCDFFile
        from GeoEco.Datasets.HDF4 import HDF4SDSCollection

        datasetsFound = []

        for url in self._QueryResults:
            queryableAttributeValues = self._GetQueryableAttributeValuesForFile(url)
            if queryableAttributeValues is None:
                continue

            if parsedExpression is not None:
                try:
                    matches = parsedExpression.eval(queryableAttributeValues)
                except Exception, e:
                    continue        # TODO: report better message
            else:
                matches = True

            if matches:
                if self._NetCDF4Available:
                    collection = NetCDFFile(url.split('/')[-1], parentCollection=self, queryableAttributeValues=queryableAttributeValues, lazyPropertyValues=self._GetLazyPropertyValuesForFile(queryableAttributeValues, self._NetCDF4Available), cacheDirectory=self.CacheDirectory)
                else:
                    collection = HDF4SDSCollection(url.split('/')[-1], parentCollection=self, queryableAttributeValues=queryableAttributeValues, lazyPropertyValues=self._GetLazyPropertyValuesForFile(queryableAttributeValues, self._NetCDF4Available), cacheDirectory=self.CacheDirectory)

                datasetsFound.extend(collection._QueryDatasets(parsedExpression, progressReporter, options, queryableAttributeValues))

        return datasetsFound

    def _GetOldestDataset(self, parsedExpression, options, parentAttrValues, dateTimeAttrName):
        return self._GetOldestOrNewestDataset(parsedExpression, options, False)

    def _GetNewestDataset(self, parsedExpression, options, parentAttrValues, dateTimeAttrName):
        return self._GetOldestOrNewestDataset(parsedExpression, options, True)

    def _GetOldestOrNewestDataset(self, parsedExpression, options, newest):

        # If we have not queried the server yet, do it now.

        if self._QueryResults is None:
            self._QueryServer()

        # Go through the list of URLs returned by the server in time
        # order, testing whether each one matches the query
        # expression. For each match, construct a NetCDFFile or
        # HDF4SDSCollection instance and query it. As soon as that
        # query returns an Grid instance, return it.

        from GeoEco.Datasets.NetCDF import NetCDFFile
        from GeoEco.Datasets.HDF4 import HDF4SDSCollection

        for url in sorted(self._QueryResults, cmp=lambda x, y: cmp(x[1:], y[1:]), reverse=newest):
            queryableAttributeValues = self._GetQueryableAttributeValuesForFile(url)
            if queryableAttributeValues is None:
                continue

            if parsedExpression is not None:
                try:
                    matches = parsedExpression.eval(queryableAttributeValues)
                except Exception, e:
                    continue        # TODO: report better message
            else:
                matches = True

            if matches:
                if self._NetCDF4Available:
                    collection = NetCDFFile(url.split('/')[-1], parentCollection=self, queryableAttributeValues=queryableAttributeValues, lazyPropertyValues=self._GetLazyPropertyValuesForFile(queryableAttributeValues, self._NetCDF4Available), cacheDirectory=self.CacheDirectory)
                else:
                    collection = HDF4SDSCollection(url.split('/')[-1], parentCollection=self, queryableAttributeValues=queryableAttributeValues, lazyPropertyValues=self._GetLazyPropertyValuesForFile(queryableAttributeValues, self._NetCDF4Available), cacheDirectory=self.CacheDirectory)

                datasets = collection._QueryDatasets(parsedExpression, None, options, queryableAttributeValues)
                if len(datasets) > 0:
                    return datasets[0]

        # We did not find a matching dataset. Return None.

        return None

    def _GetLocalFile(self, pathComponents):

        # We need a place to cache the downloaded file. Check whether
        # we or our parent collections have a cache directory defined.
        # If so, use it. If not, create a temporary one.

        cacheDirectory = None
        obj = self
        while obj is not None:
            if obj.CacheDirectory is not None:
                cacheDirectory = obj.CacheDirectory
                break
            obj = obj.ParentCollection
        
        if cacheDirectory is None:
            self.CacheDirectory = self._CreateTempDirectory()
            cacheDirectory = self.CacheDirectory

        # If the file does not already exist, download it.

        localFile = os.path.join(cacheDirectory, pathComponents[0])
        if not os.path.isfile(localFile):
            url = 'http://oceandata.sci.gsfc.nasa.gov/cgi/getfile/' + pathComponents[0]
            self._LogDebug(_(u'%(class)s 0x%(id)08X: Downloading %(url)s to %(file)s'), {u'class': self.__class__.__name__, u'id': id(self), u'url': url, u'file': localFile})

            if self._Http is None:
                from GeoEco.httplib2 import Http
                self._Http = Http(timeout=self._Timeout)
                self._Http.follow_all_redirects = True
                self._Http.force_exception_to_status_code = False
                self._RegisterForCloseAtExit()

            try:
                self._Http.download_file_with_retry(url, localFile, max_retry_time=self._MaxRetryTime, logger=self._GetLogger(), message=_(u'Failed to download file %(url)s from the NASA OceanColor server to %(file)s. The HTTP request failed with %%(e)s: %%(msg)s. Retrying...') % {u'url': url, u'file': localFile})
            except Exception, e:
                if e.__class__.__name__ == 'ExecuteAbort' and str(e).lower() == 'cancelled function':
                    raise
                if self._MaxRetryTime is not None:
                    raise RuntimeError(_(u'Failed to download file %(url)s from the NASA OceanColor server to %(file)s. The download was retried for %(retry)i seconds without success. Check that the server is operating properly, that your computer can connect to it, that you have write access to the destination directory, and that disk is not full. If necessary, contact the server\'s operator for assistance. If the server and network are operating properly, this problem could be a programming error in this tool. If you suspect one, contact the author of this tool for assistance. Error details: %(e)s: %(msg)s') % {u'url': url, u'file': localFile, u'retry': self._MaxRetryTime, u'e': e.__class__.__name__, u'msg': UserPreferredEncodingToUnicode(e)})
                raise RuntimeError(_(u'Failed to download file %(url)s from the NASA OceanColor server to %(file)s. Check that the server is operating properly, that your computer can connect to it, that you have write access to the destination directory, and that disk is not full. If necessary, contact the server\'s operator for assistance. If the server and network are operating properly, this problem could be a programming error in this tool. If you suspect one, contact the author of this tool for assistance. Error details: %(e)s: %(msg)s') % {u'url': url, u'file': localFile, u'e': e.__class__.__name__, u'msg': UserPreferredEncodingToUnicode(e)})
                
        return localFile, True          # True indicates that it is ok for the caller to delete the downloaded file after decompressing it, to save space

    def _QueryServer(self):
        # If we have not queried the server yet, do it now.

        if self._QueryResults is None:

            # Check whether the netCDF4 module is available. If it is,
            # we will access the netCDF-4 files NASA began publishing
            # in summer 2015. If not, we'll fall back to the HDF files
            # NASA published historically.

            from GeoEco.Datasets.NetCDF import NetCDFDependency

            d = NetCDFDependency()
            d.Initialize()
            if NetCDFDependency.GetNetCDFModuleName() == 'netCDF3':
                self._LogWarning(_(u'In order to access the netCDF-4 files that NASA OceanColor began publishing in June 2015, this tool requires the netCDF4 module for Python %(major)s.%(minor)s, which does not appear to be installed. Until you install it, this tool will attempt to access the HDF files that NASA used to publish, which may be out of date or no longer exist. If you are an ArcGIS user, you can upgrade to ArcGIS 10.3 or later; ArcGIS 10.3 and later automatically install the netCDF4 module when ArcGIS installs Python. If this will not work for you, you can manually install the package yourself. Go to https://pypi.python.org/pypi/netCDF4/ and look for "Windows binary installers". You must restart ArcGIS after installing it.') % {u'major': platform.python_version_tuple()[0], u'minor': platform.python_version_tuple()[1]})
                self._NetCDF4Available = False
            else:
                self._NetCDF4Available = True

            # Build the dictionary of query parameters.

            params = {'search': '',
                      'dtype': 'L3m',
                      'std_only': '1',
                      'results_as_file': '1',
                      'addurl' : '1'}

            if self._TemporalResolution is not None:
                if self._TemporalResolution.lower() == u'daily':
                    params['search'] += '*_DAY'
                elif self._TemporalResolution.lower() == u'8day':
                    params['search'] += '*_8D'
                elif self._TemporalResolution.lower() == u'monthly':
                    params['search'] += '*_MO'
                elif self._TemporalResolution.lower() == u'annual':
                    params['search'] += '*_YR'
                else:
                    raise ValueError(_(u'Programming error in this tool: The temporal resolution "%(temporalResolution)s" is not currently supported. Please contact the author of this tool for assistance.') % {u'temporalResolution': self._TemporalResolution})

            if self._Product is not None:
                params['search'] += '*_' + self._Product.lower()

            if self._NetCDF4Available:
                params['search'] += '*.nc'
            else:
                params['search'] += '*.bz2'

            if self._StartDate is not None:
                params['sdate'] = self._StartDate.strftime('%Y-%m-%d')
            else:
                params['sdate'] = ''

            if self._EndDate is not None:
                params['edate'] = self._EndDate.strftime('%Y-%m-%d')
            else:
                params['edate'] = ''

            if self._Sensor is None:
                params['sensor'] = 'all'
            elif self._Sensor.lower() == 'seawifs':
                params['sensor'] = 'seawifs'
            elif self._Sensor.lower() == 'aqua':
                params['sensor'] = 'aqua'
            elif self._Sensor.lower() == 'terra':
                params['sensor'] = 'terra'
            elif self._Sensor.lower() == 'czcs':
                params['sensor'] = 'czcs'
            elif self._Sensor.lower() == 'octs':
                params['sensor'] = 'octs'
            elif self._Sensor.lower() == 'meris':
                params['sensor'] = 'meris'
            else:
                raise ValueError(_(u'Programming error in this tool: The sensor "%(sensor)s" is not currently supported. Please contact the author of this tool for assistance.') % {u'sensor': self._Sensor})

            # POST the query to the server.

            url = 'http://oceandata.sci.gsfc.nasa.gov/search/file_search.cgi'

            self._LogDebug(_(u'%(class)s 0x%(id)08X: Querying %(url)s with parameters: %(params)s'), {u'class': self.__class__.__name__, u'id': id(self), u'url': url, u'params': repr(params)})

            if self._Http is None:
                from GeoEco.httplib2 import Http
                self._Http = Http(timeout=self._Timeout)
                self._Http.follow_all_redirects = True
                self._Http.force_exception_to_status_code = False
                self._RegisterForCloseAtExit()

            from urllib import urlencode

            try:
                response, content = self._Http.request_with_retry(url + '?' + urlencode(params), 'GET', max_retry_time=self._MaxRetryTime, logger=self._GetLogger(), message=UnicodeToUserPreferredEncoding(_(u'Failed to query the NASA OceanColor server at %(url)s. The HTTP request failed with %%(e)s: %%(msg)s. Retrying...') % {u'url': url}))
            except Exception, e:
                if e.__class__.__name__ == 'ExecuteAbort' and str(e).lower() == 'cancelled function':
                    raise
                if self._MaxRetryTime is not None:
                    raise RuntimeError(_(u'Failed to query the NASA OceanColor server at %(url)s. The query was retried for %(retry)i seconds without success. Check that the server is operating properly and that your computer can connect to it. If necessary, contact the server\'s operator for assistance. If the server and network are operating properly, this problem could be a programming error in this tool. If you suspect one, contact the author of this tool for assistance. Error details: %(e)s: %(msg)s') % {u'url': url, u'retry': self._MaxRetryTime, u'e': e.__class__.__name__, u'msg': UserPreferredEncodingToUnicode(e)})
                raise RuntimeError(_(u'Failed to query the NASA OceanColor server at %(url)s. Check that the server is operating properly and that your computer can connect to it. If necessary, contact the server\'s operator for assistance. If the server and network are operating properly, this problem could be a programming error in this tool. If you suspect one, contact the author of this tool for assistance. Error details: %(e)s: %(msg)s') % {u'url': url, u'e': e.__class__.__name__, u'msg': UserPreferredEncodingToUnicode(e)})

            # Parse the results.

            self._QueryResults = []

            lines = content.split('\n')
            for line in lines:
                if line.strip().startswith('http://oceandata.sci.gsfc.nasa.gov/cgi/getfile'):
                    self._QueryResults.append(line.strip())

            self._LogDebug(_(u'%(class)s 0x%(id)08X: %(count)s URLs were returned.'), {u'class': self.__class__.__name__, u'id': id(self), u'count': len(self._QueryResults)})

    @classmethod
    def _GetQueryableAttributeValuesForFile(cls, url):
        queryableAttributeValues = {}

        fileName = url.split('/')[-1]

        if fileName.startswith('O2_'):
            return None                     # When searching for OCTS, the server erroneously also returns files from the Ocean Colour Monitor 2. Ignore these.

        if fileName[0] in ['A', 'a']:
            queryableAttributeValues['Sensor'] = u'aqua'
        elif fileName[0] in ['T', 't']:
            queryableAttributeValues['Sensor'] = u'terra'
        elif fileName[0] in ['S', 's']:
            queryableAttributeValues['Sensor'] = u'seawifs'
        elif fileName[0] in ['C', 'C']:
            queryableAttributeValues['Sensor'] = u'czcs'
        elif fileName[0] in ['O', 'o']:
            queryableAttributeValues['Sensor'] = u'octs'
        elif fileName[0] in ['M', 'm']:
            queryableAttributeValues['Sensor'] = u'meris'
        else:
            return None

        queryableAttributeValues['SensorCode'] = fileName[0].lower()

        if '.' not in fileName or len(fileName.split('.')[1].split('_')) < 4 or fileName.split('.')[1].split('_')[0].lower() != 'l3m':
            return None

        components = fileName.split('.')[1].split('_')

        if components[1].upper() == 'DAY':
            queryableAttributeValues['TemporalResolution'] = u'daily'
        elif components[1].upper() == '8D':
            queryableAttributeValues['TemporalResolution'] = u'8day'
        elif components[1].upper() == 'MO':
            queryableAttributeValues['TemporalResolution'] = u'monthly'
        elif components[1].upper() == 'YR':
            queryableAttributeValues['TemporalResolution'] = u'annual'
        else:
            return None

        queryableAttributeValues['TemporalResolutionCode'] = components[1].upper()

        if components[-1].lower() not in ['4', '9', '4km', '9km']:
            return None

        queryableAttributeValues['SpatialResolution'] = components[-1][0] + 'km'
        queryableAttributeValues['SpatialResolutionCode'] = components[-1][0]
        
        queryableAttributeValues['ProductCode'] = u'_'.join(components[2:-1])
        
        if '_' not in queryableAttributeValues['ProductCode']:
            return None         # At the time this tool was written (December 2010, updated January 2015), all sensors used a naming convention that had at least one _ in the name. Ignore files that do not have this; they are old, or are SST, which we do not support.

        dt = datetime.datetime(int(fileName[1:5]), 1, 1) + datetime.timedelta(days=int(fileName[5:8]) - 1)
        queryableAttributeValues['DateTime'] = dt
        queryableAttributeValues['Year'] = dt.year
        queryableAttributeValues['Month'] = dt.month
        queryableAttributeValues['Day'] = dt.day
        queryableAttributeValues['Hour'] = dt.hour
        queryableAttributeValues['Minute'] = dt.minute
        queryableAttributeValues['Second'] = dt.second
        queryableAttributeValues['DayOfYear'] = int(dt.strftime('%j'))

        return queryableAttributeValues

    @classmethod
    def _GetLazyPropertyValuesForFile(cls, queryableAttributeValues, netCDF4Available):
        values = {'SpatialReference': Dataset.ConvertSpatialReference('proj4', '+proj=latlong +ellps=WGS84 +datum=WGS84 +no_defs', 'obj'),
                  'Dimensions': 'yx',
                  'Shape': {'4km': (4320, 8640), '9km': (2160, 4320)}[queryableAttributeValues['SpatialResolution']],
                  'CoordDependencies': (None, None),
                  'CoordIncrements': {'4km': (180./4320,360./8640), '9km': (180./2160,360./4320)}[queryableAttributeValues['SpatialResolution']],
                  'TIncrement': {'daily': 1, '8day': 8, 'monthly': 1, u'annual': 1}[queryableAttributeValues['TemporalResolution']],
                  'TIncrementUnit': {'daily': 'day', '8day': 'day', 'monthly': 'month', u'annual': 'year'}[queryableAttributeValues['TemporalResolution']],
                  'TSemiRegularity': {'daily': None, '8day': 'annual', 'monthly': None, u'annual': None}[queryableAttributeValues['TemporalResolution']],
                  'TCountPerSemiRegularPeriod': {'daily': None, '8day': 46, 'monthly': None, u'annual': None}[queryableAttributeValues['TemporalResolution']],
                  'TCornerCoordType': 'min',
                  'TOffsetFromParsedTime': None,
                  'CornerCoords': {'4km': (-90 + 180./4320/2, -180 + 360./8640/2), '9km': (-90 + 180./2160/2, -180 + 360./4320/2)}[queryableAttributeValues['SpatialResolution']],
                  'PhysicalDimensions': 'yx',
                  'PhysicalDimensionsFlipped': (True, False),
                  'UnscaledDataType': u'float32',
                  'UnscaledNoDataValue': -32767,
                  'ScaledDataType': None,
                  'ScaledNoDataValue': None,
                  'ScalingFunction': None,
                  'UnscalingFunction': None}

        if netCDF4Available:
            values['VariableNames'] = [queryableAttributeValues['ProductCode'].split('_',1)[-1]]

        else:
            values['SDSNames'] = ['l3m_data']
            values['SDSIndices'] = [0]

            # In 2010, Aqua, OCTS, SeaWiFS, and Terra were reprocessed and
            # the resulting HDF files all used 32-bit floats with a NoData
            # value of -32767, except for LAND_NDVI for SeaWiFS, which
            # still used uint16 with a scaling equation and a NoData value
            # of 65535.
            #
            # UPDATE October 2013: A user reported that the zeu_kpar
            # product is also still uint16 with a scaling equation, and I
            # then discovered zeu_zeul and zeu_zeum have the same problem.
            # For now, I'm going to add cases in here for these, but the
            # proper solution is to either modify this or
            # GeoEco.Datasets.HDF4.HDF4SDS to recognize the HDF-EOS
            # attributes and build a scaling equation on the fly. That
            # would be resilient to NASA adding products or changing
            # scaling equations. The current approach below is NOT
            # resilient, and will fail whenever NASA makes changes
            # (although it works just fine so long as NASA continues to
            # publish things as unscaled float32).
            #
            # UPDATE January 2015: CZCS has now been reprocessed and uses
            # 32-bit floats, like the other sensors.
            
            import numpy

            if queryableAttributeValues['ProductCode'].lower() in [u'land_ndvi', u'zeu_kpar', u'zeu_zeul', u'zeu_zeum']:
                values['UnscaledDataType'] = u'uint16'
                values['UnscaledNoDataValue'] = 65535
                values['ScaledDataType'] = u'float32'
                if queryableAttributeValues['ProductCode'].lower() == u'land_ndvi':
                    values['ScaledNoDataValue'] = numpy.cast['float32'](1.4728232599736657e-005*65535 - 0.05)
                    values['ScalingFunction'] = lambda data: numpy.cast['float32'](1.4728232599736657e-005*data - 0.05)
                    values['UnscalingFunction'] = lambda data: numpy.cast['uint16'](numpy.round((0.05 + data)/1.4728232599736657e-005))
                elif queryableAttributeValues['ProductCode'].lower() == u'zeu_kpar':
                    values['ScaledNoDataValue'] = numpy.cast['float32'](10.**(0.000043*65535 - 2.))
                    values['ScalingFunction'] = lambda data: numpy.cast['float32'](10.**(0.000043*data - 2.))
                    values['UnscalingFunction'] = lambda data: numpy.cast['uint16'](numpy.round((numpy.log10(data) + 2.)/0.000043))
                elif queryableAttributeValues['ProductCode'].lower() in [u'zeu_zeul', u'zeu_zeum']:
                    values['ScaledNoDataValue'] = numpy.cast['float32'](0.002670*65535 + 5.)
                    values['ScalingFunction'] = lambda data: numpy.cast['float32'](0.002670*data + 5.)
                    values['UnscalingFunction'] = lambda data: numpy.cast['uint16'](numpy.round((5. + data)/0.002670))
        
        return values


class OceanColorLevel3SMITimeSeries(TimeSeriesGridStack):
    __doc__ = DynamicDocString()

    def __init__(self, sensor, temporalResolution, spatialResolution, product, timeout=60, maxRetryTime=300, cacheDirectory=None):
        self.__doc__.Obj.ValidateMethodInvocation()

        # Initialize our properties.

        self._DisplayName = _(u'%(sensor)s %(temporalResolution)s %(spatialResolution)s %(product)s from NASA OceanColor') % {u'sensor': sensor, u'temporalResolution': temporalResolution, u'spatialResolution': spatialResolution, u'product': product}

        # Construct a OceanColorLevel3SMIFileSearcher for the product
        # requested by the caller.

        collection = OceanColorLevel3SMIFileSearcher(sensor, temporalResolution, product, timeout=timeout, maxRetryTime=maxRetryTime, cacheDirectory=cacheDirectory)

        # Create an expression that selects the time series of
        # datasets requested by the caller. Include all four
        # parameters as query terms. This is partially redundant
        # because we passed three of them into the
        # OceanColorLevel3SMIFileSearcher constructor above. But we do
        # it so that if there are no datasets matching the caller's
        # request, our base class will raise a descriptive exception
        # that includes the query expression.

        expression = u"Sensor = '%s' AND TemporalResolution = '%s' AND SpatialResolution = '%s' AND ProductCode = '%s'" % (sensor, temporalResolution, spatialResolution, product)

        # Initialize the base class.

        super(OceanColorLevel3SMITimeSeries, self).__init__(collection, expression=expression, reportProgress=False)

    def _GetDisplayName(self):
        return self._DisplayName

    @classmethod
    def _GetTimeCoordsFromQueryableAttributeValues(cls, queryableAttributeValues):
        try:
            if u'DateTime' in queryableAttributeValues and isinstance(queryableAttributeValues[u'DateTime'], datetime.datetime):
                startDate = queryableAttributeValues[u'DateTime']
                temporalResolution = None
                
                if u'TemporalResolution' in queryableAttributeValues and isinstance(queryableAttributeValues[u'TemporalResolution'], basestring) and queryableAttributeValues[u'TemporalResolution'] in [u'daily', '8day', 'monthly', 'annual']:
                    temporalResolution = queryableAttributeValues[u'TemporalResolution']
                elif u'TemporalResolutionCode' in queryableAttributeValues and isinstance(queryableAttributeValues[u'TemporalResolutionCode'], basestring) and queryableAttributeValues[u'TemporalResolutionCode'] in [u'DAY', '8D', 'MO', 'YR']:
                    temporalResolution = queryableAttributeValues[u'TemporalResolutionCode']
                else:
                    return [startDate, None, None]
                    
                if temporalResolution in [u'daily', u'DAY']:
                    return [startDate, startDate + datetime.timedelta(0.5), startDate + datetime.timedelta(1) - datetime.timedelta(seconds=1)]
                
                if temporalResolution in [u'8day', u'8D']:
                    if startDate.month == 12 and startDate.day >= 24:
                        endDate = datetime.datetime(startDate.year + 1, 1, 1)
                        return [startDate, startDate + (endDate - startDate) / 2, endDate - datetime.timedelta(seconds=1)]
                    else:
                        return [startDate, startDate + datetime.timedelta(4), startDate + datetime.timedelta(8) - datetime.timedelta(seconds=1)]
                    
                if temporalResolution in [u'monthly', u'MO']:
                    if startDate.month == 12:
                        endDate = datetime.datetime(startDate.year + 1, 1, 1)
                    else:
                        endDate = datetime.datetime(startDate.year, startDate.month + 1, 1)
                    return [startDate, startDate + (endDate - startDate) / 2, endDate - datetime.timedelta(seconds=1)]

                if temporalResolution in [u'annual', u'YR']:
                    return [startDate, datetime.datetime(startDate.year, 7, 1), datetime.datetime(startDate.year + 1, 1, 1) - datetime.timedelta(seconds=1)]

            return [None, None, None]
        except:
            return [None, None, None]

    @classmethod
    def CreateArcGISRasters(cls, sensor, temporalResolution, spatialResolution, product,
                            outputWorkspace, mode=u'add', rasterNameExpressions=[u'%(Sensor)s', u'%(TemporalResolution)s', u'%(SpatialResolution)s', u'%(ProductCode)s', u'%%Y', u'%(SensorCode).1s%%Y%%j%(EndDate)s.L3m_%(TemporalResolutionCode).3s_%(ProductCode)s_%(SpatialResolution)s.img'], rasterCatalog=None, cacheDirectory=None,
                            rotationOffset=None, spatialExtent=None, startDate=None, endDate=None, 
                            timeout=60, maxRetryTime=300,
                            calculateStatistics=True, buildPyramids=False):
        cls.__doc__.Obj.ValidateMethodInvocation()
        grid = cls(sensor, temporalResolution, spatialResolution, product, timeout, maxRetryTime, cacheDirectory)
        try:
            from GeoEco.Datasets.ArcGIS import ArcGISWorkspace, ArcGISRaster
            grid = cls._RotateAndClip(grid, rotationOffset, spatialExtent, startDate, endDate)
            workspace = ArcGISWorkspace(outputWorkspace, ArcGISRaster, pathCreationExpressions=rasterNameExpressions, cacheTree=True, queryableAttributes=tuple(grid.GetAllQueryableAttributes() + [QueryableAttribute(u'DateTime', _(u'Date'), DateTimeTypeMetadata())]))
            workspace.ImportDatasets(GridSliceCollection(grid).QueryDatasets(), mode, calculateStatistics=calculateStatistics, buildPyramids=buildPyramids)
            if rasterCatalog is not None:
                workspace.ToRasterCatalog(rasterCatalog, grid.GetSpatialReference(u'ArcGIS'), tQACoordType=u'min', tCoordFunction=cls._GetTimeCoordsFromQueryableAttributeValues, overwriteExisting=True)
        finally:
            grid.Close()
        return outputWorkspace

    @classmethod
    def CreateClimatologicalArcGISRasters(cls, sensor, temporalResolution, spatialResolution, product,
                                          statistic, binType,
                                          outputWorkspace, mode=u'add', rasterNameExpressions=[u'%(Sensor)s', u'%(TemporalResolution)s', u'%(SpatialResolution)s', u'%(ProductCode)s', u'%(ClimatologyBinType)s_Climatology', u'%(Sensor)s_%(ProductCode)s_%(ClimatologyBinName)s_%(Statistic)s.img'], cacheDirectory=None,
                                          binDuration=1, startDayOfYear=1,
                                          rotationOffset=None, spatialExtent=None, startDate=None, endDate=None,
                                          timeout=60, maxRetryTime=300,
                                          calculateStatistics=True, buildPyramids=False):
        cls.__doc__.Obj.ValidateMethodInvocation()
        grid = cls(sensor, temporalResolution, spatialResolution, product, timeout, maxRetryTime, cacheDirectory)
        try:
            from GeoEco.Datasets.ArcGIS import ArcGISWorkspace, ArcGISRaster
            grid = cls._RotateAndClip(grid, rotationOffset, spatialExtent, startDate, endDate)
            collection = ClimatologicalGridCollection(grid, statistic, binType, binDuration, startDayOfYear, reportProgress=True)
            workspace = ArcGISWorkspace(outputWorkspace, ArcGISRaster, pathCreationExpressions=rasterNameExpressions, cacheTree=True, queryableAttributes=tuple(collection.GetAllQueryableAttributes()))
            workspace.ImportDatasets(collection.QueryDatasets(), mode, calculateStatistics=calculateStatistics, buildPyramids=buildPyramids)
        finally:
            grid.Close()
        return outputWorkspace

    @classmethod
    def InterpolateAtArcGISPoints(cls, sensor, temporalResolution, spatialResolution, product,
                                   points, valueField, tField, method=u'Nearest', cacheDirectory=None, where=None, noDataValue=None,
                                   timeout=60, maxRetryTime=300,
                                   orderByFields=None, numBlocksToCacheInMemory=256, xBlockSize=16, yBlockSize=16, tBlockSize=3):
        cls.__doc__.Obj.ValidateMethodInvocation()
        grid = cls(sensor, temporalResolution, spatialResolution, product, timeout, maxRetryTime, cacheDirectory)
        try:
            if orderByFields is not None:
                orderBy = u', '.join(map(lambda f: f + u' ASC', orderByFields))
            else:
                from GeoEco.ArcGIS import GeoprocessorManager
                if GeoprocessorManager.GetArcGISMajorVersion() > 9 or GeoprocessorManager.GetArcGISMinorVersion() >= 2:
                    orderBy = tField + u' ASC'
                else:
                    orderBy = None
            from GeoEco.Datasets.ArcGIS import ArcGISTable
            from GeoEco.SpatialAnalysis.Interpolation import Interpolator
            Interpolator.InterpolateGridsValuesForTableOfPoints([grid], ArcGISTable(points), [valueField], tField=tField, where=where, orderBy=orderBy, method=method, noDataValue=noDataValue, gridsWrap=True, numBlocksToCacheInMemory=numBlocksToCacheInMemory, xBlockSize=xBlockSize, yBlockSize=yBlockSize, tBlockSize=tBlockSize)
        finally:
            grid.Close()
        return points

    @classmethod
    def _RotateAndClip(cls, grid, rotationOffset, spatialExtent, startDate, endDate):
        if rotationOffset is not None:
            grid = RotatedGlobalGrid(grid, rotationOffset, u'Map units')

        xMin, yMin, xMax, yMax = None, None, None, None
        if spatialExtent is not None:
            from GeoEco.Types import EnvelopeTypeMetadata
            xMin, yMin, xMax, yMax = EnvelopeTypeMetadata.ParseFromArcGISString(spatialExtent)

        if spatialExtent is not None or startDate is not None or endDate is not None:
            if startDate is not None:
                startDate = datetime.datetime(startDate.year, startDate.month, startDate.day, 0, 0, 0)
            if endDate is not None:
                endDate = datetime.datetime(endDate.year, endDate.month, endDate.day, 23, 59, 59)
                
            grid = ClippedGrid(grid, u'Map coordinates', xMin=xMin, xMax=xMax, yMin=yMin, yMax=yMax, tMin=startDate, tMax=endDate)
            
        return grid


class OceanColorLevel2Converter(object):
    __doc__ = DynamicDocString()

    @classmethod
    def L2LacSstHdfToArcGISPoints(cls, hdf, points, variable=u'sst', spatialExtent=None,
                                  qualityLevel=1, recoverBadPixels=False, checkATMFAIL=True, checkLAND=True, checkPRODWARN=False, checkHIGLINT=False, checkHILT=True, checkHISATZEN=True, checkCOASTZ=False, checkSTRAYLIGHT=True, checkCLDICE=True, checkCOCCOLITH=False, checkTURBIDW=False, checkHISOLZEN=True, checkLOWLW=True, checkCHLFAIL=False, checkNAVWARN=True, checkABSAER=False, checkMAXAERITER=True, checkMODGLINT=False, checkCHLWARN=False, checkATMWARN=True, checkSEAICE=True, checkNAVFAIL=True, checkFILTER=False, checkHIPOL=False, checkPRODFAIL=False,
                                  overwriteExisting=False):
        cls.__doc__.Obj.ValidateMethodInvocation()
        [lon, lat, data] = cls._ReadL2SSTHDF(hdf, variable, qualityLevel, recoverBadPixels, checkATMFAIL, checkLAND, checkPRODWARN, checkHIGLINT, checkHILT, checkHISATZEN, checkCOASTZ, checkSTRAYLIGHT, checkCLDICE, checkCOCCOLITH, checkTURBIDW, checkHISOLZEN, checkLOWLW, checkCHLFAIL, checkNAVWARN, checkABSAER, checkMAXAERITER, checkMODGLINT, checkCHLWARN, checkATMWARN, checkSEAICE, checkNAVFAIL, checkFILTER, checkHIPOL, checkPRODFAIL, spatialExtent)
        cls._NumpyArraysToArcGISPoints(points, variable, lon, lat, data)

    @classmethod
    def L2LacSstHdfToArcGISRaster(cls, hdf, raster, variable=u'sst', spatialExtent=None,
                                  qualityLevel=1, recoverBadPixels=False, checkATMFAIL=True, checkLAND=True, checkPRODWARN=False, checkHIGLINT=False, checkHILT=True, checkHISATZEN=True, checkCOASTZ=False, checkSTRAYLIGHT=True, checkCLDICE=True, checkCOCCOLITH=False, checkTURBIDW=False, checkHISOLZEN=True, checkLOWLW=True, checkCHLFAIL=False, checkNAVWARN=True, checkABSAER=False, checkMAXAERITER=True, checkMODGLINT=False, checkCHLWARN=False, checkATMWARN=True, checkSEAICE=True, checkNAVFAIL=True, checkFILTER=False, checkHIPOL=False, checkPRODFAIL=False,
                                  templateRaster=None, power=1., searchRadius=u'FIXED', distance=1154.7, numPoints=0, barrier=None,
                                  overwriteExisting=False):

        # Read the HDF into numpy arrays.

        [lon, lat, data] = cls._ReadL2SSTHDF(hdf, variable, qualityLevel, recoverBadPixels, checkATMFAIL, checkLAND, checkPRODWARN, checkHIGLINT, checkHILT, checkHISATZEN, checkCOASTZ, checkSTRAYLIGHT, checkCLDICE, checkCOCCOLITH, checkTURBIDW, checkHISOLZEN, checkLOWLW, checkCHLFAIL, checkNAVWARN, checkABSAER, checkMAXAERITER, checkMODGLINT, checkCHLWARN, checkATMWARN, checkSEAICE, checkNAVFAIL, checkFILTER, checkHIPOL, checkPRODFAIL, spatialExtent)

        # Create a temporary point feature class from the HDF. If Arc
        # 10 or later, use the in-memory workspace.

        from GeoEco.ArcGIS import GeoprocessorManager
        gp = GeoprocessorManager.GetWrappedGeoprocessor()
        
        if GeoprocessorManager.GetArcGISMajorVersion() >= 10:
            tempPoints = 'in_memory\\MGET_OCL2_Points'
        else:
            from GeoEco.DataManagement.Directories import TemporaryDirectory
            tempDir = TemporaryDirectory()
            tempPoints = os.path.join(tempDir.Path, u'points.shp')

        try:
            cls._NumpyArraysToArcGISPoints(tempPoints, variable, lon, lat, data)

            # Now determine the coordinate system, cell size, and
            # extent we should use for the output raster. If the
            # caller provided a template raster, look these up from
            # it. If not, create an Albers coordinate system suitable
            # for the points, use 1000 m as the cell size, and let the
            # IDW tool determine the spatial extent.

            if templateRaster is not None:
                d = gp.Describe(templateRaster)
                outputCS = d.SpatialReference
                cellSize = d.MeanCellWidth
                [left, bottom, right, top] = EnvelopeTypeMetadata.ParseFromArcGISString(d.Extent)
            else:
                centralMeridian = (float(min(lon)) + float(max(lon))) / 2
                latRange = float(max(lat)) - float(min(lat))
                parallel1 = float(max(lat)) - latRange/6
                parallel2 = float(min(lat)) + latRange/6
                outputCS = "PROJCS['Albers_Auto_Generated',GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Albers'],PARAMETER['False_Easting',0.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',%r],PARAMETER['Standard_Parallel_1',%r],PARAMETER['Standard_Parallel_2',%r],PARAMETER['Latitude_Of_Origin',%r],UNIT['Meter',1.0]] # GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]]" % (centralMeridian, parallel1, parallel2, centralMeridian)
                cellSize = 1000.
                [left, bottom, right, top] = [None, None, None, None]

            # Configure the ArcGIS environment settings to produce a
            # raster conforming to that configuration.

            oldOutputCoordinateSystem = gp.OutputCoordinateSystem
            try:
                gp.OutputCoordinateSystem = outputCS
                oldCellSize = gp.CellSize
                try:
                    gp.CellSize = cellSize
                    oldExtent = gp.Extent
                    try:
                        if templateRaster is not None:
                            oldSnapRaster = -1
                            gp.Extent = '%r %r %r %r' % (left, bottom, right, top)
                            if GeoprocessorManager.GetArcGISMajorVersion() > 9 or GeoprocessorManager.GetArcGISMajorVersion() == 9 and GeoprocessorManager.GetArcGISMinorVersion() >= 3:
                                oldSnapRaster = gp.SnapRaster
                                gp.SnapRaster = templateRaster

                        # If teh caller wants SST, use Spatial
                        # Analyst's IDW to interpolate the raster.

                        Logger.Info(_(u'Interpolating %(out)s from %(in)s.') % {u'in': tempPoints, u'out': raster})

                        if variable == u'sst':
                            if searchRadius.upper() == 'VARIABLE':
                                idwSearchRadius = 'VARIABLE %r %r' % (numPoints, distance)
                            else:
                                idwSearchRadius = 'FIXED %r %r' % (distance, numPoints)
                            gp.Idw_sa(tempPoints, variable, raster, cellSize, power, idwSearchRadius, barrier)

                        # Otherwise--the caller wants sst_qual or
                        # l2_flags, which are integer data that should
                        # not be interpolated with Spatial Analyst's
                        # interpolators--use Point To Raster.
                        
                        else:
                            gp.PointToRaster_conversion (tempPoints, variable, raster, 'MAXIMUM', None, cellSize)

                    finally:
                        try:
                            if templateRaster is not None:
                                gp.Extent = oldExtent
                                if oldSnapRaster != -1:
                                    gp.SnapRaster = oldSnapRaster
                        except:
                            pass
                finally:
                    try:
                        gp.CellSize = oldCellSize
                    except:
                        pass
            finally:
                try:
                    gp.OutputCoordinateSystem = oldOutputCoordinateSystem
                except:
                    pass
        finally:
            try:
                if GeoprocessorManager.GetArcGISMajorVersion() >= 10 and gp.Exists(tempPoints):
                    gp.Delete_management(tempPoints)
            except:
                pass

    @classmethod
    def MultipleL2LacSstHdfsToArcGISPoints(cls, hdfs, points, variable=u'sst', spatialExtent=None, startDate=None, numDays=None,
                                           qualityLevel=1, recoverBadPixels=False, checkATMFAIL=True, checkLAND=True, checkPRODWARN=False, checkHIGLINT=False, checkHILT=True, checkHISATZEN=True, checkCOASTZ=False, checkSTRAYLIGHT=True, checkCLDICE=True, checkCOCCOLITH=False, checkTURBIDW=False, checkHISOLZEN=True, checkLOWLW=True, checkCHLFAIL=False, checkNAVWARN=True, checkABSAER=False, checkMAXAERITER=True, checkMODGLINT=False, checkCHLWARN=False, checkATMWARN=True, checkSEAICE=True, checkNAVFAIL=True, checkFILTER=False, checkHIPOL=False, checkPRODFAIL=False,
                                           overwriteExisting=False):
        #cls.__doc__.Obj.ValidateMethodInvocation()

        # For each HDF, extract the date from the file name. If the
        # date is within the range specified by startDate and numDays,
        # read the data from it.

        import numpy

        mergedLon = numpy.array([], 'float32')
        mergedLat = numpy.array([], 'float32')
        if variable.startswith('qual'):
            mergedData = numpy.array([], 'int8')
        else:
            mergedData = numpy.array([], 'float32')

        for hdf in hdfs:
            hdfDate = datetime.datetime(*(time.strptime(os.path.basename(hdf)[1:8], '%Y%j')[0:6]))
            if startDate is None or hdfDate >= startDate and hdfDate < startDate + datetime.timedelta(days=numDays):
                [lon, lat, data] = cls._ReadL2SSTHDF(hdf, variable, qualityLevel, recoverBadPixels, checkATMFAIL, checkLAND, checkPRODWARN, checkHIGLINT, checkHILT, checkHISATZEN, checkCOASTZ, checkSTRAYLIGHT, checkCLDICE, checkCOCCOLITH, checkTURBIDW, checkHISOLZEN, checkLOWLW, checkCHLFAIL, checkNAVWARN, checkABSAER, checkMAXAERITER, checkMODGLINT, checkCHLWARN, checkATMWARN, checkSEAICE, checkNAVFAIL, checkFILTER, checkHIPOL, checkPRODFAIL, spatialExtent)
                mergedLon = numpy.hstack([mergedLon, lon])
                mergedLat = numpy.hstack([mergedLat, lat])
                mergedData = numpy.hstack([mergedData, data])

        # Write all of the retrieved points to the feature class.
        
        cls._NumpyArraysToArcGISPoints(points, variable, mergedLon, mergedLat, mergedData)

    @classmethod
    def _ReadL2SSTHDF(cls, hdf, variable, qualityLevel, recoverBadPixels, checkATMFAIL, checkLAND, checkPRODWARN, checkHIGLINT, checkHILT, checkHISATZEN, checkCOASTZ, checkSTRAYLIGHT, checkCLDICE, checkCOCCOLITH, checkTURBIDW, checkHISOLZEN, checkLOWLW, checkCHLFAIL, checkNAVWARN, checkABSAER, checkMAXAERITER, checkMODGLINT, checkCHLWARN, checkATMWARN, checkSEAICE, checkNAVFAIL, checkFILTER, checkHIPOL, checkPRODFAIL, spatialExtent):

        # Read the SDS from the HDF into parallel 1D numpy arrays,
        # dropping any records for which any SDS contains NoData.

        Logger.Info(_(u'Reading L2_LAC_SST HDF file %(hdf)s.') % {u'hdf': hdf})

        if qualityLevel < 4:
            if recoverBadPixels:
                [data, qual, l2Flags, lon, lat] = cls._ReadL2HDF(hdf, [variable, 'qual_sst', 'l2_flags'])
            else:
                [data, qual, lon, lat] = cls._ReadL2HDF(hdf, [variable, 'qual_sst'])
        else:
            [data, lon, lat] = cls._ReadL2HDF(hdf, [variable])

        # If the caller is requesting the l2_flags variable, force the
        # most significant bit to zero.
        # http://oceancolor.gsfc.nasa.gov/VALIDATION/flags.html lists
        # this bit as "spare" but NASA seems to be setting it. (I
        # recall a recent document showing that they are setting
        # values for the spare bits for experimental purposes.) The
        # problem is that if this bit is set but no others are, then
        # the corresponding 32-bit signed integer is -2147483648. This
        # is the value that ArcGIS uses for the No Data value for
        # 32-bit signed integer rasters. Therefore, pixels that have
        # perfect quality (by virtue of having no l2_flags set) end up
        # being converted to No Data by L2LacSstHdfToArcGISRaster.
        # That does not make sense. They should be converted to 0. To
        # deal with this problem, just clear the most significant bit.
        # This will force all values to positive, and cause "perfect"
        # pixels to be zero.

        import numpy

        if variable == 'l2_flags':
            data = numpy.bitwise_and(data, int(2**31-1))

        # Apply the minimum qualityLevel, if requested.

        # TODO: check for whether it is daytime or nighttime, and whether it is Aqua

        if qualityLevel < 4:
            if recoverBadPixels:
                useFlags = 0
                useFlags += checkATMFAIL * 2**(1-1)
                useFlags += checkLAND * 2**(2-1)
                useFlags += checkPRODWARN * 2**(3-1)
                useFlags += checkHIGLINT * 2**(4-1)
                useFlags += checkHILT * 2**(5-1)
                useFlags += checkHISATZEN * 2**(6-1)
                useFlags += checkCOASTZ * 2**(7-1)
                useFlags += checkSTRAYLIGHT * 2**(9-1)
                useFlags += checkCLDICE * 2**(10-1)
                useFlags += checkCOCCOLITH * 2**(11-1)
                useFlags += checkTURBIDW * 2**(12-1)
                useFlags += checkHISOLZEN * 2**(13-1)
                useFlags += checkLOWLW * 2**(15-1)
                useFlags += checkCHLFAIL * 2**(16-1)
                useFlags += checkNAVWARN * 2**(17-1)
                useFlags += checkABSAER * 2**(18-1)
                useFlags += checkMAXAERITER * 2**(20-1)
                useFlags += checkMODGLINT * 2**(21-1)
                useFlags += checkCHLWARN * 2**(22-1)
                useFlags += checkATMWARN * 2**(23-1)
                useFlags += checkSEAICE * 2**(25-1)
                useFlags += checkNAVFAIL * 2**(26-1)
                useFlags += checkFILTER * 2**(27-1)
                useFlags += checkHIPOL * 2**(30-1)
                useFlags += checkPRODFAIL * 2**(31-1)
                isValid = (qual <= qualityLevel) | ((qual <= 3) & (numpy.bitwise_and(l2Flags, useFlags) == 0))
            else:
                isValid = qual <= qualityLevel
            data = data[isValid]
            lon = lon[isValid]
            lat = lat[isValid]

        # Clip to the spatial extent, if requested.

        if spatialExtent is not None:
            from GeoEco.Types import EnvelopeTypeMetadata
            lonMin, latMin, lonMax, latMax = EnvelopeTypeMetadata.ParseFromArcGISString(spatialExtent)
            
            isInRegion = (lon >= lonMin) & (lon <= lonMax) & (lat >= latMin) & (lat <= latMax)
            data = data[isInRegion]
            lon = lon[isInRegion]
            lat = lat[isInRegion]

        # Return successfully.

        return [lon, lat, data]

    @classmethod
    def _ReadL2HDF(cls, hdfPath, sdsNames):
        import numpy

        # Open the HDF file and validate that we can process the
        # requested variables.

        from GeoEco.Datasets.HDF4 import HDF4SDSCollection

        hdf = HDF4SDSCollection(hdfPath, lazyPropertyValues={'Dimensions': u'yx', 'PhysicalDimensions': u'yx', 'PhysicalDimensionsFlipped': (False, False)})

        if len(hdf.QueryDatasets('SDSName = \'sst4\'', reportProgress=False)) > 0:      # Automatically change the variables from sst and qual_sst to sst4 and qual_sst4, if necessary.
            newSdsNames = []
            for i in range(len(sdsNames)):
                if sdsNames[i] == 'sst':
                    newSdsNames.append('sst4')
                elif sdsNames[i] == 'qual_sst':
                    sdsNames[i] = 'qual_sst4'
                    newSdsNames.append('qual_sst4')
                else:
                    newSdsNames.append(sdsNames[i])
            sdsNames = newSdsNames

        if len(hdf.QueryDatasets('SDSName = \'Longitude\'', reportProgress=False)) > 0:
            sdsNames.append('Longitude')
        else:
            sdsNames.append('longitude')

        if len(hdf.QueryDatasets('SDSName = \'Latitude\'', reportProgress=False)) > 0:
            sdsNames.append('Latitude')
        else:
            sdsNames.append('latitude')

        variableSDSes = []
        lonIndex = None
        latIndex = None
        dataIndex = None
        for i in range(len(sdsNames)):
            sdsName = sdsNames[i]
            datasets = hdf.QueryDatasets('SDSName = \'%s\'' % sdsName, reportProgress=False)
            
            if len(datasets) <= 0:
                raise ValueError(_('The HDF file %(hdf)s does not contain a variable named %(sdsName)s. Please check the file to see if it is compatible with this tool. If you have any doubts, contact the MGET development team.') % {u'hdf': hdfPath, u'sdsName': sdsName})
            if len(datasets[0].Shape) != 2:
                raise ValueError(_('The variable %(sdsName)s in HDF file %(hdf)s has %(dim)i dimensions. It must have 2 dimensions to be processed by this tool. Please check the file to see if it is compatible with this tool. If you have any doubts, contact the MGET development team.') % {u'hdf': hdfPath, u'sdsName': sdsName, u'dim': datasets[0].Shape})
            
            variableSDSes.append(datasets[0])

            if sdsName.lower() == 'longitude':
                lonIndex = i
            elif sdsName.lower() == 'latitude':
                latIndex = i
            elif dataIndex is None:
                dataIndex = i
    
        # For each variable, look up the No Data values and scaling
        # equations and apply them.

        for sds in variableSDSes:
            sdsObj = sds._GetSDS()

            # First try to look up the No Data value that should be
            # used before any scaling equation is applied. Note that
            # NASA GSFC OceanColor uses bad_value_scaled (not
            # bad_value_unscaled) for this value. See
            # http://oceancolor.gsfc.nasa.gov/REPROCESSING/R2009/format/'
            # for more information.

            for attrName in ['_FillValue', 'bad_value_scaled']:
                if hasattr(sdsObj, attrName):
                    if sds.DataType.startswith('float'):
                        sds.SetLazyPropertyValue('UnscaledNoDataValue', int(getattr(sdsObj, attrName)))
                    else:
                        sds.SetLazyPropertyValue('UnscaledNoDataValue', float(getattr(sdsObj, attrName)))
                    break

            # Now, if we found that the variable contains integers,
            # try to look up the scaling equation for converting those
            # integers to floats. (When we say "scaled" here, we mean
            # the float value computed from the integers, not the
            # other way around, as NASA GSFC OceanColor seems to mean
            # in their documentation.)

            if not sds.DataType.startswith('float'):
                for (slopeAttrName, interceptAttrName) in {'slope': 'intercept', 'scale_factor': 'add_offset'}.items():
                    if hasattr(sdsObj, slopeAttrName) and hasattr(sdsObj, interceptAttrName):
                        slope = float(getattr(sdsObj, slopeAttrName))
                        intercept = float(getattr(sdsObj, interceptAttrName))
                        if slope != 1. or intercept != 0:
                            sds.SetLazyPropertyValue('ScaledDataType', u'float32')
                            def lambdaFactory(slope2, intercept2):                              # See http://stackoverflow.com/questions/938429/scope-of-python-lambda-functions-and-their-parameters for an explanation of this
                                return lambda d: numpy.cast['float32'](slope2*d + intercept2)
                            sds.SetLazyPropertyValue('ScalingFunction', lambdaFactory(slope, intercept))
                            if sds.UnscaledNoDataValue is not None:
                                sds.SetLazyPropertyValue('ScaledNoDataValue', sds.GetLazyPropertyValue('ScalingFunction')(sds.UnscaledNoDataValue))
                            break

        # Read the data for each variable.

        data = []
        for sds in variableSDSes:
            data.append(sds.Data[:])

        # Check the latitude and longitude arrays. If they have fewer
        # rows than the data variable, it means that the coordinate
        # variables were subsampled to save space in the HDF file, as
        # described in this forum post:
        # http://oceancolor.gsfc.nasa.gov/forum/oceancolor/topic_show.pl?tid=758
        # If that happens, we have to interpolate the missing values
        # using a cubic spline interpolator.

        cntl_pt_cols = None
        tempDir = None

        for i in [lonIndex, latIndex]:
            if data[i].shape[0] != data[dataIndex].shape[0]:
                raise ValueError(_('In HDF file %(hdf)s, the %(coord)s variable does not have the same number of rows as the %(data)s variable. The two variables have dimensions %(dim1)r and %(dim2)r, respectively. This is unexpected. The file cannot be processed by this tool. Please contact the MGET development team for assistance') % {u'hdf': hdfPath, u'coord': sdsNames[i], u'data': sdsNames[dataIndex], u'dim1': data[i].shape, u'dim2': data[dataIndex].shape})

            if data[i].shape[1] != data[dataIndex].shape[1]:
                Logger.Info(_(u'Interpolating missing rows of the %(coord)s variable. (NASA subsampled this variable to save space in the HDF file.)') % {u'coord': sdsNames[i]})

                # The number of rows are too small. For example, there
                # might be 170 rows instead of the expected 1354. Read
                # the cntl_pt_cols variable. This will be a scalar
                # vector with 170 values, which specify the 1-based
                # indices (e.g. between 1 and 1354) that each of the
                # 170 rows corresponds to. cntl_pt_cols[0] will always
                # be 1; cntl_pt_cols[-1] will always be the highest
                # possible index (1354). cntl_pt_cols typically has a
                # stride of 8, according to the NASA forum post. We
                # need to interpolate rows for the missing indices.
                #
                # Use MATLAB to do the interpolation. Before settling
                # on this approach I considered a few pure Python
                # spline interplators, a Cython based interpolator,
                # and scipy. The pure Python interpolators were too
                # slow. The Cython interpolator would require adding
                # build support for Cython; too much work. Scipy would
                # require taking a dependency on it. Since MGET
                # already requires the MATLAB MCR for several things,
                # I thought it was better to use that existing
                # dependency rather than take a new one. (I also
                # considered extracting my own copy of
                # scipy.interpolate for MGET's private use, but it
                # looked too complicated to build for this problem,
                # requiring a compilation of linpack from fortran
                # files, or something like that.)
                #
                # We prefer to call MATLAB directly right here but
                # there is a continuing incompatibility between MATLAB
                # DLLs and ArcGIS DLLs (they both try to load their
                # own incompatible versions of xerces-c_2_7.dll) so we
                # have to do it in a separate process.

                if cntl_pt_cols is None:
                    hdf._Open()
                    cntl_pt_cols = hdf._HDF.select('cntl_pt_cols')[:]
                    
                    from GeoEco.DataManagement.Directories import TemporaryDirectory
                    tempDir = TemporaryDirectory()

                    cntl_pt_cols_file = os.path.join(tempDir.Path, 'cntl_pt_cols.dat')
                    numpy.cast['float32'](cntl_pt_cols).tofile(cntl_pt_cols_file)

                coordFile = os.path.join(tempDir.Path, sdsNames[i] + '.dat')
                data[i].tofile(coordFile)
                
                outputFile = os.path.join(tempDir.Path, sdsNames[i] + '_interpolated.dat')

                from GeoEco.DataManagement.Processes import ChildProcess
                import win32api

                args = [os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), 'SpatialAnalysis', 'InterpolateSpline1D.py'),
                        cntl_pt_cols_file,
                        coordFile,
                        repr(data[i].shape[0]),
                        repr(data[i].shape[1]),
                        outputFile,
                        str(data[i].dtype.name)]

                oldLogInfoAsDebug = Logger.GetLogInfoAsDebug()
                Logger.SetLogInfoAsDebug(True)
                try:
                    ChildProcess.ExecuteProgram(ChildProcess.GetPythonExecutable(),
                                                arguments=args,
                                                windowState=u'invisible',
                                                maxRunTime=None)
                finally:
                    Logger.SetLogInfoAsDebug(oldLogInfoAsDebug)

                data[i] = numpy.fromfile(outputFile, data[i].dtype.name).reshape(data[dataIndex].shape)

        # Flatten the arrays for each variable.

        for i in range(len(data)):
            data[i] = data[i].flatten()

        # If any variables have No Data values, remove rows for which
        # any variable is No Data.

        if any([sds.NoDataValue is not None for sds in variableSDSes]):
            hasData = numpy.array([True] * len(data[0]))
            
            for i in range(len(data)):
                if variableSDSes[i].NoDataValue is not None:
                    hasData = numpy.logical_or(hasData, data[i] == variableSDSes[i].NoDataValue)

            for i in range(len(data)):
                data[i] = data[i][hasData]

        # Return successfully.

        return data

    @classmethod
    def _NumpyArraysToArcGISPoints(cls, points, variable, lon, lat, data):

        # Create the output point feature class.

        from GeoEco.Datasets.ArcGIS import ArcGISWorkspace, ArcGISTable

        workspace = ArcGISWorkspace(os.path.dirname(points), ArcGISTable, pathCreationExpressions=['%(TableName)s'], queryableAttributes=(QueryableAttribute('TableName', 'Table name', UnicodeStringTypeMetadata()),))
        table = workspace.CreateTable(os.path.basename(points), geometryType='Point', spatialReference=Dataset.ConvertSpatialReference('proj4', '+proj=latlong +ellps=WGS84 +datum=WGS84 +no_defs', 'obj'))
        if data.dtype.name in ['int8', 'uint8']:
            table.AddField(variable, 'int16')
        else:
            table.AddField(variable, data.dtype.name)

        # Insert the points.

        Logger.Info(_(u'Inserting %(rowCount)i points into %(points)s.') % {u'rowCount': len(data), 'points': table.DisplayName})

        isFloat = data.dtype.name.startswith('float')
        ogr = ArcGISTable._ogr()
        cur = table.OpenInsertCursor(rowCount=len(data))
        try:
            for i in range(len(data)):
                point = ogr.Geometry(ogr.wkbPoint)
                point.SetPoint_2D(0, float(lon[i]), float(lat[i]))
                cur.SetGeometry(point)
                if isFloat:
                    cur.SetValue(variable, float(data[i]))
                else:
                    cur.SetValue(variable, int(data[i]))
                cur.InsertRow()
        finally:
            del cur


###############################################################################
# Metadata: module
###############################################################################

from GeoEco.ArcGIS import ArcGISDependency, ArcGISExtensionDependency
from GeoEco.Dependencies import PythonAggregatedModuleDependency
from GeoEco.Datasets.ArcGIS import _CalculateStatisticsDescription, _BuildPyramidsDescription
from GeoEco.Metadata import *
from GeoEco.Types import *

AddModuleMetadata(shortDescription=_(u'DatasetCollections and Grids that represent products published by the NASA GSFC OceanColor Group.'))

###############################################################################
# Metadata: OceanColorLevel3SMIFileSearcher class
###############################################################################

_OceanColorLevel3SMI_LongDescription = _(
u"""The `NASA Goddard Space Flight Center (GSFC) OceanColor Group <http://oceancolor.gsfc.nasa.gov/>`_
publishes a variety of satellite image products derived from ocean
color observations made by polar-orbiting sensors such as MODIS,
SeaWiFS, OCTS, and CZCS. The most popular product is an estimate of
chlorophyll-a concentration.

This %(name)s accesses the Level 3 Standard Mapped Image (SMI)
products, which have global spatial extent, use a geographic
coordinate system with the WGS 1984 datum, and have square cells with
either 1/12 or 1/24 degree resolution (about 9.3 km or 4.6 km at the
equator).

NASA publishes the SMI products as collections of compressed HDF
version 4 files that are downloadable from the OceanColor web site.
This %(name)s automatically downloads, decompresses, and reads HDF
files as they are needed. Unless you specify a directory to cache the
files, they will be stored in your user TEMP directory and deleted
when processing is finished.

**References**

To cite the use of NASA OceanColor data in a publication, please see
`these instructions <http://oceancolor.gsfc.nasa.gov/forum/oceancolor/topic_show.pl?tid=474>`_.

For a list of publications from the NASA OceanColor Group, see
`this page <http://oceancolor.gsfc.nasa.gov/cgi/obpgpubs.cgi>`_.""")

AddClassMetadata(OceanColorLevel3SMIFileSearcher,
    shortDescription=_(u'A DatasetCollection that queries NASA GSFC for Level 3 Standard Mapped Images (SMI).'),
    longDescription=_OceanColorLevel3SMI_LongDescription % {u'name': 'class'})

# Public constructor: OceanColorLevel3SMIFileSearcher.__init__

AddMethodMetadata(OceanColorLevel3SMIFileSearcher.__init__,
    shortDescription=_(u'OceanColorLevel3SMIFileSearcher constructor.'),
    isExposedToPythonCallers=True)

AddArgumentMetadata(OceanColorLevel3SMIFileSearcher.__init__, u'self',
    typeMetadata=ClassInstanceTypeMetadata(cls=OceanColorLevel3SMIFileSearcher),
    description=_(u'%s instance.') % OceanColorLevel3SMIFileSearcher.__name__)

AddArgumentMetadata(OceanColorLevel3SMIFileSearcher.__init__, u'sensor',
    typeMetadata=UnicodeStringTypeMetadata(canBeNone=True, allowedValues=[u'Aqua', u'CZCS', u'MERIS', u'OCTS', u'SeaWiFS', u'Terra'], makeLowercase=True),
    description=_(
u"""Sensor to use, one of:

* Aqua - the Moderate Resolution Imaging Spectroradiometer (MODIS)
  sensor carried by the Aqua satellite. Aqua datasets start in July
  2002 and were still being collected at the time this tool was
  written.

* CZCS - the Coastal Zone Color Scanner (CZCS) carried by the Nimbus 7
  satellite. CZCS datasets start in September 1978 and end in June
  1986.

* MERIS - the Medium Resolution Imaging Spectrometer (MERIS) carried
  by the ENVISAT-1 satellite. MERIS datasets start in May 2010 and end
  in April 2012.

* OCTS - the Ocean Color and Temperature Scanner (OCTS) carried by the
  ADEOS-1 satellite. OCTS datasets start in November 1996 and end in
  June 1997. Although the mission was designed to last several years,
  ADEOS-1 stopped communicating after nine months due to the failure
  of its solar power system.

* SeaWiFS - Sea-viewing Wide Field-of-view Sensor (SeaWiFS) carried by
  the SeaStar satellite. SeaWiFS datasets start in September 1997 and
  were still being collected at the time this tool was written.
  SeaWiFS has been operating far beyond its designed lifetime and has
  experienced periodic failures in recent years. In particular, as of
  this writing, little or no data are available for the time periods
  of January to March 2008, July and August 2008, late April to mid
  July 2009, and September to November 2009.

* Terra - the Moderate Resolution Imaging Spectroradiometer (MODIS)
  sensor carried by the Terra satellite. Terra datasets start in
  February 2000 and were still being collected at the time this tool
  was written. **Warning:** Due to problems with the sensor scan
  mirror, ocean color observations from MODIS Terra are considered to
  be significantly less accurate than those from MODIS Aqua or
  SeaWiFS, and NASA recommends
  `here <http://oceancolor.gsfc.nasa.gov/forum/oceancolor/topic_show.pl?tid=3734>`_
  that "if you have a choice between any other sensor and MODIS Terra,
  choose the other sensor." NASA devised statistical algorithms to
  correct the data somewhat; for more information, see Franz et al.
  (2008) and Kwiatkowska et al. (2008).

The NASA OceanColor Group may publish data for other sensors, but they
are not supported by this tool at this time. If you need one of those
products, please contact the author of this tool to see if support may
be added.

**References**

Kwiatkowska, E.J., B.A. Franz, G. Meister, C.R. McClain, and X. Xiong
(2008). Cross-Calibration of ocean color bands from Moderate
Resolution Imaging Spectroradiometer on Terra platform. Applied Optics
47(36): 6796-6810.

Franz, B.A., E.J. Kwiatkowska, G. Meister, and C.R. McClain (2008).
Moderate Resolution Imaging Spectroradiometer on Terra: limitations
for ocean color applications, Journal of Applied Remote Sensing 2:
023525.
"""),
    arcGISDisplayName=_(u'Sensor'))

##* VIIRS - the Visible Infrared Imaging Radiometer Suite (VIIRS)
##  carried by the Suomi National Polar-orbiting Partnership (Suomi NPP)
##  satellite, previously known as the National Polar-orbiting
##  Operational Environmental Satellite System Preparatory Project (NPP)
##  or NPP-Bridge. VIIRS datasets start in Januaray 2012 and were still
##  being collected at the time this tool was written.


AddArgumentMetadata(OceanColorLevel3SMIFileSearcher.__init__, u'temporalResolution',
    typeMetadata=UnicodeStringTypeMetadata(canBeNone=True, allowedValues=[u'Daily', u'8day', u'Monthly', u'Annual'], makeLowercase=True),
    description=_(
u"""Temporal resolution to use, one of:

* Daily - daily images. There are 365 during normal years and 366
  during leap years.

* 8day - 8-day images. There are 46 per year. The first image of the
  year starts on January 1. The duration of the last image of the year
  is five days during normal years and six days during leap years.

* Monthly - monthly images.

* Annual - annual images.

Although NASA may publish MODIS SST images at other temporal
resolutions, they are not supported at this time. If you need one of
those products, please contact the author of this tool to see if
support may be added.

The ocean color sensors experience occasional transient failures that
prevent data from being collected, sometimes for an extended period.
NASA opted not to produce any images for these periods. These missing
images are represented as time slices filled with the NoData value.
For example, during 2004, NASA produced only 43 8-day images of
chlorophyll-a concentration for the Aqua satellite. Thus, of the 46
8-day time slices for 2004, 43 have some valid pixels while 3 are
filled entirely with the NoData value."""),
    arcGISDisplayName=_(u'Temporal resolution'))

AddArgumentMetadata(OceanColorLevel3SMIFileSearcher.__init__, u'product',
    typeMetadata=UnicodeStringTypeMetadata(canBeNone=True, makeLowercase=False),    # October 2015: Product names are case sensitive: they appear as variable names in the new netCDF4 files
    description=_(
u"""Product code of the NASA Level 3 Standard Mapped Image (SMI)
product to use, such as CHL_chlor_a for chlorophyll concentration.

The product code is case sensitive. For example, for chlorophyll you
must specify CHL_chlor_a, not chl_chlor_a, CHL_CHLOR_A, and so on.

The products that are available depend on the sensor. Newer sensors
such as SeaWiFS and MODIS provide more products. The product must be
specified using a code assigned by NASA. Most users will be interested
in the chlorophyll-a concentration product, which has the code
CHL_chlor_a for all sensors.

For all sensors, NASA provides a set of "standard" products that are
well tested and believed to be of wide interest. For a few sensors,
NASA also provides "evaluation" and "test" products, which are less
well-tested and of narrower interest. Please see NASA documentation
for more information on the products you are interested in.

Here, we list all of the products we were aware of when this tool was
developed. If you are aware of product that is not listed here, you
may try its product code. The product code is defined by the
characters that appear in NASA's file name between the temporal
resolution and spatial resolution codes. For example, in the file
O1997164.L3m_DAY_CHL_chlor_a_9km.bz2, the product code is
CHL_chlor_a.

This tool only supports L3 SMI products at 4 km and 9 km resolution.
It does not support L0, L1, or L2 products. For those, please try the
`SeaDAS tool <http://seadas.gsfc.nasa.gov/>`_. It does not support
binned products, or products at other spatial resolutions.

**Aqua and Terra MODIS - Standard Products:**

Most MODIS products are available at both 9 km and 4 km resolution.

* CDOM_cdom_index - Chromorphic dissolved organic matter index
* CHL_chlor_a - Chlorophyll-a concentration (mg m-3)
* FLH_ipar - Instantaneous photosynthetically available radiation (Einstein / m2 / sec)
* FLH_nflh - Normalized flourescence line height (mW / cm2 / um / sr)
* KD490_Kd_490 - Diffuse attenuation coefficient at 490 nm (m-1)
* PAR_par - Photosynthetically available radiation (Einstein / m2 / day)
* PIC_pic - Particulate inorganic carbon (mol / m3)
* POC_poc - Particulate organic carbon (mol / m3)
* RRS_Rrs_412 - Remote sensing reflectance at 412 nm (sr-1)
* RRS_Rrs_443 - Remote sensing reflectance at 443 nm (sr-1)
* RRS_Rrs_469 - Remote sensing reflectance at 469 nm (sr-1)
* RRS_Rrs_488 - Remote sensing reflectance at 488 nm (sr-1)
* RRS_Rrs_531 - Remote sensing reflectance at 531 nm (sr-1)
* RRS_Rrs_547 - Remote sensing reflectance at 547 nm (sr-1)
* RRS_Rrs_555 - Remote sensing reflectance at 555 nm (sr-1)
* RRS_Rrs_645 - Remote sensing reflectance at 645 nm (sr-1)
* RRS_Rrs_667 - Remote sensing reflectance at 667 nm (sr-1)
* RRS_Rrs_678 - Remote sensing reflectance at 678 nm (sr-1)
* RRS_angstrom - Angstrom coefficient
* RRS_aot_869 - Aerosol optical thickness at 869 nm

**Aqua MODIS - Evaluation Products:**

* GSM_adg_443_gsm - Absorption due to gelbstof and detritus at 443 nm (GSM) (m-1)
* GSM_bbp_443_gsm - Particulate backscatter at 443 nm (GSM) (m-1)
* GSM_chl_gsm - Chlorophyll-a concentration (GSM) (mg m-3)
* KDLEE_Kd_412_lee - Diffuse attenuation at 412 nm (Lee) (m-1)
* KDLEE_Kd_443_lee - Diffuse attenuation at 443 nm (Lee) (m-1)
* KDLEE_Kd_488_lee - Diffuse attenuation at 488 nm (Lee) (m-1)
* KDLEE_Zeu_lee - Euphotic depth (Lee) (m)
* QAA_a_443_qaa - Total absorption at 443 nm (QAA) (m-1)
* QAA_adg_443_qaa - Absorption due to gelbstof and detritus at 443 nm (QAA) (m-1)
* QAA_aph_443_qaa - Absorption due to phytoplankton at 443 nm (QAA) (m-1)
* QAA_bbp_443_qaa - Particulate backscatter at 443 nm (QAA) (m-1)
* ZEU_KPAR - Diffuse attenuation coefficient for PAR (KPAR, Morel) (m-1)
* ZEU_ZEUL - Euphotic depth (Lee) (m)
* ZEU_ZEUM - Euphotic depth (Morel) (m)

**Aqua MODIS - Test Products:**

* GIOP01_a_443_giop - Total absorption at 443 nm (m-1)
* GIOP01_a_547_giop - Total absorption at 547 nm (m-1)
* GIOP01_adg_443_giop - Absorption due to gelbstof and detritus at 443 nm (m-1)
* GIOP01_adg_s_giop - Spectral slope for gelbstof and detrital absorption
* GIOP01_aph_443_giop - Absorption due to phytoplankton at 443 nm (m-1)
* GIOP01_aph_547_giop - Absorption due to phytoplankton at 547 nm (m-1)
* GIOP01_bb_443_giop - Total backscatter at 443 nm (m-1)
* GIOP01_bb_547_giop - Total backscatter at 547 nm (m-1)
* GIOP01_bbp_443_giop - Particulate backscatter at 443 nm (m-1)
* GIOP01_bbp_s_giop - Spectral slope for particulate backscatter
* GIOP01_chl_giop - Chlorophyll-a concentration (mg m-3)
* GIOP01_rrsdiff_giop - Relative remote sensing reflectance difference

**CZCS - Standard Products:**

CZCS products are available at both 9 km and 4 km resolution.

* CHL_chlor_a - Chlorophyll-a concentration (mg m-3)
* KD490_Kd_490 - Diffuse attenuation coefficient at 490 nm (m-1)
* RRS_aot_670 - Aerosol optical thickness at 670 nm
* RRS_Rrs_443 - Remote sensing reflectance at 443 nm (sr-1)
* RRS_Rrs_520 - Remote sensing reflectance at 520 nm (sr-1)
* RRS_Rrs_550 - Remote sensing reflectance at 550 nm (sr-1)
* RRS_Rrs_670 - Remote sensing reflectance at 670 nm (sr-1)

**MERIS - Standard Products:**

MERIS products are available at both 9 km and 4 km resolution.

* CHL_chlor_a - Chlorophyll-a concentration (mg m-3)
* KD490_Kd_490 - Diffuse attenuation coefficient at 490 nm (m-1)
* PAR_par - Photosynthetically available radiation (Einstein / m2 / day)
* PIC_pic - Particulate inorganic carbon (mol / m3)
* POC_poc - Particulate organic carbon (mol / m3)
* RRS_angstrom - Angstrom coefficient
* RRS_aot_865 - Aerosol optical thickness at 865 nm
* RRS_Rrs_413 - Remote sensing reflectance at 413 nm (sr-1)
* RRS_Rrs_443 - Remote sensing reflectance at 443 nm (sr-1)
* RRS_Rrs_490 - Remote sensing reflectance at 490 nm (sr-1)
* RRS_Rrs_510 - Remote sensing reflectance at 510 nm (sr-1)
* RRS_Rrs_560 - Remote sensing reflectance at 560 nm (sr-1)
* RRS_Rrs_620 - Remote sensing reflectance at 620 nm (sr-1)
* RRS_Rrs_665 - Remote sensing reflectance at 665 nm (sr-1)
* RRS_Rrs_681 - Remote sensing reflectance at 681 nm (sr-1)
* RRS_Rrs_709 - Remote sensing reflectance at 709 nm (sr-1)

**MERIS - Evaluation Products:**

There are many of these available; we are not going to enumerate them
here. Please see the NASA OceanColor Level 3 Browser website for a
list of products.

**OCTS - Standard Products:**

OCTS products are available at only at 9 km resolution.

* CHL_chlor_a - Chlorophyll-a concentration (mg m-3)
* KD490_Kd_490 - Diffuse attenuation coefficient at 490 nm (m-1)
* PIC_pic - Particulate inorganic carbon (mol / m3)
* RRS_angstrom - Angstrom coefficient
* RRS_aot_862 - Aerosol optical thickness at 862 nm
* RRS_Rrs_412 - Remote sensing reflectance at 412 nm (sr-1)
* RRS_Rrs_443 - Remote sensing reflectance at 443 nm (sr-1)
* RRS_Rrs_490 - Remote sensing reflectance at 490 nm (sr-1)
* RRS_Rrs_516 - Remote sensing reflectance at 416 nm (sr-1)
* RRS_Rrs_565 - Remote sensing reflectance at 565 nm (sr-1)
* RRS_Rrs_667 - Remote sensing reflectance at 667 nm (sr-1)

**SeaWiFS - Standard Products:**

All SeaWiFS products are available only at 9 km resolution, except for
LAND_NDVI, which is also available at 4 km.

* CDOM_cdom_index - Chromorphic dissolved organic matter index
* CHL_chlor_a - Chlorophyll-a concentration (mg m-3)
* KD490_Kd_490 - Diffuse attenuation coefficient at 490 nm (m-1)
* LAND_NDVI - Normalized difference vegetation index
* PAR_par - Photosynthetically available radiation (Einstein / m2 / day)
* PIC_pic - Particulate inorganic carbon (mol / m3)
* POC_poc - Particulate organic carbon (mol / m3)
* RRS_angstrom - Angstrom coefficient
* RRS_aot_865 - Aerosol optical thickness at 865 nm
* RRS_Rrs_412 - Remote sensing reflectance at 412 nm (sr-1)
* RRS_Rrs_443 - Remote sensing reflectance at 443 nm (sr-1)
* RRS_Rrs_490 - Remote sensing reflectance at 490 nm (sr-1)
* RRS_Rrs_510 - Remote sensing reflectance at 410 nm (sr-1)
* RRS_Rrs_555 - Remote sensing reflectance at 555 nm (sr-1)
* RRS_Rrs_670 - Remote sensing reflectance at 670 nm (sr-1)

**SeaWiFS - Evaluation Products:**

* GSM_adg_443_gsm - Absorption due to gelbstof and detritus at 443 nm (GSM) (m-1)
* GSM_bbp_443_gsm - Particulate backscatter at 443 nm (GSM) (m-1)
* GSM_chl_gsm - Chlorophyll-a concentration (GSM) (mg m-3)
* KDLEE_Kd_412_lee - Diffuse attenuation at 412 nm (Lee) (m-1)
* KDLEE_Kd_443_lee - Diffuse attenuation at 443 nm (Lee) (m-1)
* KDLEE_Kd_490_lee - Diffuse attenuation at 490 nm (Lee) (m-1)
* KDLEE_Zeu_lee - Euphotic depth (Lee) (m)
* QAA_a_443_qaa - Total absorption at 443 nm (QAA) (m-1)
* QAA_adg_443_qaa - Absorption due to gelbstof and detritus at 443 nm (QAA) (m-1)
* QAA_aph_443_qaa - Absorption due to phytoplankton at 443 nm (QAA) (m-1)
* QAA_bbp_443_qaa - Particulate backscatter at 443 nm (QAA) (m-1)

**SeaWiFS - Test Products:**

* GIOP01_a_443_giop - Total absorption at 443 nm (m-1)
* GIOP01_a_555_giop - Total absorption at 555 nm (m-1)
* GIOP01_adg_443_giop - Absorption due to gelbstof and detritus at 443 nm (m-1)
* GIOP01_adg_s_giop - Spectral slope for gelbstof and detrital absorption
* GIOP01_aph_443_giop - Absorption due to phytoplankton at 443 nm (m-1)
* GIOP01_aph_555_giop - Absorption due to phytoplankton at 555 nm (m-1)
* GIOP01_bb_443_giop - Total backscatter at 443 nm (m-1)
* GIOP01_bb_555_giop - Total backscatter at 555 nm (m-1)
* GIOP01_bbp_443_giop - Particulate backscatter at 443 nm (m-1)
* GIOP01_bbp_s_giop - Spectral slope for particulate backscatter
* GIOP01_chl_giop - Chlorophyll-a concentration (mg m-3)
* GIOP01_rrsdiff_giop - Relative remote sensing reflectance difference
"""),
    arcGISDisplayName=_(u'Level 3 SMI product code'))

AddArgumentMetadata(OceanColorLevel3SMIFileSearcher.__init__, u'startDate',
    typeMetadata=DateTimeTypeMetadata(canBeNone=True),
    description=_(
u"""Start date of the date range to search. The time component of the
start date is ignored and assumed to be 00:00:00.

The file search operation is implemented by the NASA OceanColor
server. The server searches files according to the first day of the
time period they apply to. For example, if you search for files
between 2-January-2005 and 4-January-2005, only files with daily
temporal resolution will be returned, because none of the 8-day,
monthly, or annual files begin on those dates. On the other hand, if
you search for files between 1-January-2005 and 4-January-2005, then
8-day, monthly, and annual files will also be returned because they
all begin on January 1."""))

AddArgumentMetadata(OceanColorLevel3SMIFileSearcher.__init__, u'endDate',
    typeMetadata=DateTimeTypeMetadata(canBeNone=True),
    description=_(
u"""End date of the date range to search. The time component of the
end date is ignored and assumed to be 23:59:59."""))

AddArgumentMetadata(OceanColorLevel3SMIFileSearcher.__init__, u'timeout',
    typeMetadata=IntegerTypeMetadata(minValue=1, canBeNone=True),
    description=_(
u"""Number of seconds to wait for the server to respond before failing
with a timeout error.

If you also provide a Maximum Retry Time and it is larger than the
timeout value, the failed request will be retried automatically (with
the same timout value) until it succeeds or the Maximum Retry Time has
elapsed.

If you receive a timeout error you should investigate the server to
determine if it is malfunctioning or just slow. Check the OceanColor
website to see if NASA has posted a notice about the problem, or
contact the NASA directly. If the server just slow, increase the
timeout value to a larger number, to give the server more time to
respond."""),
    arcGISDisplayName=_(u'Timeout value'),
    arcGISCategory=_(u'Network options'))

AddArgumentMetadata(OceanColorLevel3SMIFileSearcher.__init__, u'maxRetryTime',
    typeMetadata=IntegerTypeMetadata(minValue=1, canBeNone=True),
    description=_(
u"""Number of seconds to retry requests to the server before giving
up.

Use this parameter to cope with transient failures. For example, you
may find that the server is rebooted nightly during a maintenance
cycle. If you start a long running operation and want it to run
overnight without failing, set the maximum retry time to a duration
that is longer than the time that the server is offline during the
maintenance cycle.

To maximize performance while minimizing load during failure
situations, retries are scheduled with progressive delays:

* The first retry is issued immediately.

* Then, so long as fewer than 10 seconds have elapsed since the
  original request was issued, retries are issued every second.

* After that, retries are issued every 30 seconds until the maximum
  retry time is reached or the request succeeds.
"""),
    arcGISDisplayName=_(u'Maximum retry time'),
    arcGISCategory=_(u'Network options'))

AddArgumentMetadata(OceanColorLevel3SMIFileSearcher.__init__, u'cacheDirectory',
    typeMetadata=DirectoryTypeMetadata(canBeNone=True),
    description=_(
u"""Directory for caching local copies of downloaded files. A cache
directory is optional but highly recommended if you plan to repeatedly
access data for the same range of dates.

NASA partitions ocean color data into collections of compressed HDF
files according to the sensor, temporal resolution, spatial
resolution, product code, and date. These files have global spatial
extent and typically range from 5 to 70 MB in size. Thus, even if you
are only interested in a small region of the planet--even just a
single point location--this tool must still download a global file
each time slice that is needed. This can take a long time if many
files are needed.

When this tool needs a file, it will first check the cache directory
to see if the file was downloaded and cached during a prior run. If it
was, data will be read directly from that file. If not, the file will
be downloaded, decompressed, and stored in the cache directory for
later use.

If you use a cache directory, be aware of these common pitfalls:

* The caching algorithm permits the directory to grow to infinite size
  and never deletes any cached files. If you access a large number of
  files (e.g. 10 years of daily images) they will all be added to the
  cache. Be careful that you do not fill up your hard disk. To
  mitigate this, manually delete the entire cache or specific files
  within it when they are no longer needed.

* The caching algorithm stores uncompressed files so that they may be
  accessed quickly, without incuring a decompression step every time
  they are needed. To save space on your hard disk, we highly
  recommend you enable compression of the cache directory by the
  operating system. In Windows Explorer, right click on the directory,
  select Properties, click Advanced, and enable "Compress contents to
  save disk space".

* Due to limitations in the caching algorithm, it cannot detect when
  NASA reprocesses data products and replaces files on the server with
  updated versions, thereby making the cached files obsolete. Thus, if
  NASA republishes a product with improved data values, the caching
  algorithm will continue to use the old, obsolete values. To mitigate
  this, you should monitor when NASA reprocesses their products and
  delete the cached files when they become obsolete.
"""),
    arcGISDisplayName=_(u'Cache directory'))

AddResultMetadata(OceanColorLevel3SMIFileSearcher.__init__, u'collection',
    typeMetadata=ClassInstanceTypeMetadata(cls=OceanColorLevel3SMIFileSearcher),
    description=_(u'%s instance.') % OceanColorLevel3SMIFileSearcher.__name__)

###############################################################################
# Metadata: OceanColorLevel3SMITimeSeries class
###############################################################################

AddClassMetadata(OceanColorLevel3SMITimeSeries,
    shortDescription=_(u'Time series of Level 3 Standard Mapped Images (SMI) published by the NASA GSFC OceanColor Group.'),
    longDescription=_OceanColorLevel3SMI_LongDescription % {u'name': 'class'})

# Public constructor: OceanColorLevel3SMITimeSeries.__init__

AddMethodMetadata(OceanColorLevel3SMITimeSeries.__init__,
    shortDescription=_(u'OceanColorLevel3SMITimeSeries constructor.'),
    isExposedToPythonCallers=True)

AddArgumentMetadata(OceanColorLevel3SMITimeSeries.__init__, u'self',
    typeMetadata=ClassInstanceTypeMetadata(cls=OceanColorLevel3SMITimeSeries),
    description=_(u'%s instance.') % OceanColorLevel3SMITimeSeries.__name__)

CopyArgumentMetadata(OceanColorLevel3SMIFileSearcher.__init__, u'sensor', OceanColorLevel3SMITimeSeries.__init__, u'sensor')
OceanColorLevel3SMITimeSeries.__init__.__doc__.Obj.GetArgumentByName(u'sensor').Type.CanBeNone = False

CopyArgumentMetadata(OceanColorLevel3SMIFileSearcher.__init__, u'temporalResolution', OceanColorLevel3SMITimeSeries.__init__, u'temporalResolution')
OceanColorLevel3SMITimeSeries.__init__.__doc__.Obj.GetArgumentByName(u'temporalResolution').Type.CanBeNone = False

AddArgumentMetadata(OceanColorLevel3SMITimeSeries.__init__, u'spatialResolution',
    typeMetadata=UnicodeStringTypeMetadata(allowedValues=[u'4km', u'9km'], makeLowercase=True),
    description=_(
u"""Spatial resolution to use, one of:

* 4km - the grid has a cell size of 1/24 geographic degree, or about
  4.64 km at the equator, with 8640 columns and 4320 rows. 

* 9km - the grid has a cell size of 1/12 geographic degree, or about
  9.28 km at the equator, with 4320 columns and 2160 rows.
"""),
    arcGISDisplayName=_(u'Spatial resolution'))

CopyArgumentMetadata(OceanColorLevel3SMIFileSearcher.__init__, u'product', OceanColorLevel3SMITimeSeries.__init__, u'product')
OceanColorLevel3SMITimeSeries.__init__.__doc__.Obj.GetArgumentByName(u'product').Type.CanBeNone = False

CopyArgumentMetadata(OceanColorLevel3SMIFileSearcher.__init__, u'timeout', OceanColorLevel3SMITimeSeries.__init__, u'timeout')
CopyArgumentMetadata(OceanColorLevel3SMIFileSearcher.__init__, u'maxRetryTime', OceanColorLevel3SMITimeSeries.__init__, u'maxRetryTime')
CopyArgumentMetadata(OceanColorLevel3SMIFileSearcher.__init__, u'cacheDirectory', OceanColorLevel3SMITimeSeries.__init__, u'cacheDirectory')

# Public method: OceanColorLevel3SMITimeSeries.CreateArcGISRasters

AddMethodMetadata(OceanColorLevel3SMITimeSeries.CreateArcGISRasters,
    shortDescription=_(u'Creates rasters for a Level 3 Standard Mapped Image (SMI) product published by the NASA GSFC OceanColor Group.'),
    longDescription=_OceanColorLevel3SMI_LongDescription % {u'name': 'tool'},
    isExposedToPythonCallers=True,
    isExposedByCOM=True,
    isExposedAsArcGISTool=True,
    arcGISDisplayName=_(u'Create Rasters for NASA OceanColor L3 SMI Product'),
    arcGISToolCategory=_(u'Data Products\\NASA GSFC OceanColor Group\\L3 Products'),
    dependencies=[ArcGISDependency(9, 1), PythonAggregatedModuleDependency('numpy')])

AddArgumentMetadata(OceanColorLevel3SMITimeSeries.CreateArcGISRasters, u'cls',
    typeMetadata=ClassOrClassInstanceTypeMetadata(cls=OceanColorLevel3SMITimeSeries),
    description=_(u'OceanColorLevel3SMITimeSeries class or instance.'))

CopyArgumentMetadata(OceanColorLevel3SMITimeSeries.__init__, u'sensor', OceanColorLevel3SMITimeSeries.CreateArcGISRasters, u'sensor')
CopyArgumentMetadata(OceanColorLevel3SMITimeSeries.__init__, u'temporalResolution', OceanColorLevel3SMITimeSeries.CreateArcGISRasters, u'temporalResolution')
CopyArgumentMetadata(OceanColorLevel3SMITimeSeries.__init__, u'spatialResolution', OceanColorLevel3SMITimeSeries.CreateArcGISRasters, u'spatialResolution')
CopyArgumentMetadata(OceanColorLevel3SMITimeSeries.__init__, u'product', OceanColorLevel3SMITimeSeries.CreateArcGISRasters, u'product')

AddArgumentMetadata(OceanColorLevel3SMITimeSeries.CreateArcGISRasters, u'outputWorkspace',
    typeMetadata=ArcGISWorkspaceTypeMetadata(createParentDirectories=True),
    description=_(
u"""Directory or geodatabase to receive the rasters.

Unless you have a specific reason to store the rasters in a
geodatabase, we recommend you store them in a directory because it
will be much faster and allows the rasters to be organized in a tree.
If you do store the rasters in a geodatabase, you must change the
Raster Name Expressions parameter; see below for more
information."""),
    arcGISDisplayName=_(u'Output workspace'))

AddArgumentMetadata(OceanColorLevel3SMITimeSeries.CreateArcGISRasters, u'mode',
    typeMetadata=UnicodeStringTypeMetadata(allowedValues=[u'Add', u'Replace'], makeLowercase=True),
    description=_(
u"""Overwrite mode, one of:

* Add - create rasters that do not exist and skip those that already
  exist. This is the default.

* Replace - create rasters that do not exist and overwrite those that
  already exist.

'Add' should be appropriate for nearly all situations. One situation
in which 'Replace' is appropriate is when NASA reprocesses the entire
dataset with improved algorithms and releases updated images. In that
case, you may wish to recreate all of your rasters using the updated
data.

The ArcGIS Overwrite Outputs geoprocessing setting has no effect on
this tool. If 'Replace' is selected the rasters will be overwritten,
regardless of the ArcGIS Overwrite Outputs setting."""),
    arcGISDisplayName=_(u'Overwrite mode'))

AddArgumentMetadata(OceanColorLevel3SMITimeSeries.CreateArcGISRasters, u'rasterNameExpressions',
    typeMetadata=ListTypeMetadata(elementType=UnicodeStringTypeMetadata(), minLength=1),
    description=_(
u"""List of expressions specifying how the output rasters should be
named.

The default expression assumes you are storing rasters in a file
system directory and creates them in a tree structure with names that
imitate those used by NASA. When storing rasters in a directory, the
final expression specifies the file name of the raster and any
preceding expressions specify subdirectories. The extension of the
final expression determines the output raster format: .asc for ArcInfo
ASCII Grid, .bmp for BMP, .gif for GIF, .img for an ERDAS IMAGINE
file, .jpg for JPEG, .jp2 for JPEG 2000, .png for PNG, .tif for
GeoTIFF, or no extension for ArcInfo Binary Grid. The default
expression uses .img.

When storing rasters in a geodatabase, you should provide only one
expression. That expression specifies the raster's name.

Each expression may contain any sequence of characters permitted by
the output workspace. Each expression may optionally contain one or
more of the following case-sensitive codes. The tool replaces the
codes with appropriate values when creating each raster:

* %(Sensor)s - name of the sensor, either "aqua", "czcs", "octs",
  "meris", "seawifs", or "terra".

* %(SensorCode)s - abbreviation for the sensor, either "A", "C", "O",
  "S", or "T", used in NASA's file naming scheme.

* %(TemporalResolution)s - temporal resolution, either "daily",
  "8day", "monthly", or "annual".

* %(TemporalResolutionCode)s - abbreviation for the temporal
  resolution, either "DAY", "8D", "MO", or "YR", used in NASA's file
  naming scheme.

* %(SpatialResolution)s - spatial resolution, either "4km" or "9km".

* %(SpatialResolutionCode)s - abbreviation for the spatial
  resolution, either "4" or "9", used in NASA's file naming scheme for
  certain products.
  
* %(ProductCode)s - product code for the Level 3 Standard Mapped Image
  (SMI) product represented in the output raster.

* %%Y - four-digit year of the raster.

* %%m - two-digit month of the first day of the raster.

* %%d - two-digit day of the month of the first day of the raster.

* %%j - three-digit day of the year of the first day of the raster.

* %(EndDate)s - date of the last day of the raster in the format
  "YYYYjjj" where YYYY is the four-digit year and jjj is the
  three-digit day of the year. This date is used in NASA's file naming
  scheme for temporal resolutions other than "daily".
"""),
    arcGISDisplayName=_(u'Raster name expressions'))

AddArgumentMetadata(OceanColorLevel3SMITimeSeries.CreateArcGISRasters, u'rasterCatalog',
    typeMetadata=ArcGISRasterCatalogTypeMetadata(canBeNone=True, mustBeDifferentThanArguments=[u'outputWorkspace'], createParentDirectories=True),
    description=_(
u"""Raster catalog to create.

This parameter requires ArcGIS 9.3 or later.

If this parameter is specified, after the tool finishes creating
rasters, it will create an unmanaged raster catalog and import all of
the rasters in the output workspace into it. The catalog will have
fields for all of the codes specified in the Raster Name Expressions
as well as fields for the start date, center date, and end date of
each raster. You can then use the catalog to create animations using
the ArcGIS 10 Time Slider.

WARNING: The raster catalog will be deleted and re-created each time
the tool is executed. Beware of this if you plan to add your own
fields to the catalog after it has been created."""),
    direction=u'Output',
    arcGISDisplayName=_(u'Output raster catalog'),
    dependencies=[ArcGISDependency(9, 3)])

CopyArgumentMetadata(OceanColorLevel3SMITimeSeries.__init__, u'cacheDirectory', OceanColorLevel3SMITimeSeries.CreateArcGISRasters, u'cacheDirectory')

AddArgumentMetadata(OceanColorLevel3SMITimeSeries.CreateArcGISRasters, u'rotationOffset',
    typeMetadata=FloatTypeMetadata(canBeNone=True),
    description=_(
u"""Degrees to rotate the output rasters about the polar axis. If not
provided, the output rasters will be centered on the Prime Meridian
and have x coordinates ranging from -180 to 180.

Use this parameter to shift the center longitude to a different
location. Positive values shift it to the east, negative values to the
west. For example, to center the output rasters on the Pacific ocean
and use x coordinates ranging from 0 to 360 rather than -180 to 180,
provide 180 for this parameter.

The rasters can only be rotated in whole grid cells. The value you
provide will be rounded off to the closest cell."""),
    arcGISDisplayName=_(u'Rotate raster by'),
    arcGISCategory=_(u'Spatiotemporal extent'))

AddArgumentMetadata(OceanColorLevel3SMITimeSeries.CreateArcGISRasters, u'spatialExtent',
    typeMetadata=EnvelopeTypeMetadata(canBeNone=True),
    description=_(
u"""Spatial extent of the output rasters, in degrees. If not provided,
the spatial extent will be the entire planet.

This parameter is applied after the rotation parameter and uses
coordinates that result after rotation. For example, if the rotation
parameter was 180, the resulting rasters will have x coordinates
ranging from 0 to 360. The spatial extent should be expressed in those
coordinates.

The rasters can only be clipped in whole grid cells. The values you
provide will be rounded off to the closest cell."""),
    arcGISDisplayName=_(u'Spatial extent'),
    arcGISCategory=_(u'Spatiotemporal extent'))

AddArgumentMetadata(OceanColorLevel3SMITimeSeries.CreateArcGISRasters, u'startDate',
    typeMetadata=DateTimeTypeMetadata(canBeNone=True),
    description=_(
u"""Start date for the rasters to create. Rasters will be created for
images that occur on or after the start date and on or before the end
date. If the start date is not provided, the date of the oldest image
will be used.

The time component of the start date is ignored."""),
    arcGISDisplayName=_(u'Start date'),
    arcGISCategory=_(u'Spatiotemporal extent'))

AddArgumentMetadata(OceanColorLevel3SMITimeSeries.CreateArcGISRasters, u'endDate',
    typeMetadata=DateTimeTypeMetadata(canBeNone=True),
    description=_(
u"""End date for the rasters to create. Rasters will be created for
images that occur on or after the start date and on or before the end
date. If the end date is not provided, the date of the most recent
image will be used.

The time component of the end date is ignored."""),
    arcGISDisplayName=_(u'End date'),
    arcGISCategory=_(u'Spatiotemporal extent'))

CopyArgumentMetadata(OceanColorLevel3SMITimeSeries.__init__, u'timeout', OceanColorLevel3SMITimeSeries.CreateArcGISRasters, u'timeout')
CopyArgumentMetadata(OceanColorLevel3SMITimeSeries.__init__, u'maxRetryTime', OceanColorLevel3SMITimeSeries.CreateArcGISRasters, u'maxRetryTime')

AddArgumentMetadata(OceanColorLevel3SMITimeSeries.CreateArcGISRasters, u'calculateStatistics',
    typeMetadata=BooleanTypeMetadata(),
    description=_CalculateStatisticsDescription,
    arcGISDisplayName=_(u'Calculate statistics'),
    arcGISCategory=_(u'Additional raster processing options'))

AddArgumentMetadata(OceanColorLevel3SMITimeSeries.CreateArcGISRasters, u'buildPyramids',
    typeMetadata=BooleanTypeMetadata(),
    description=_BuildPyramidsDescription,
    arcGISDisplayName=_(u'Build pyramids'),
    arcGISCategory=_(u'Additional raster processing options'))

AddResultMetadata(OceanColorLevel3SMITimeSeries.CreateArcGISRasters, u'updatedOutputWorkspace',
    typeMetadata=ArcGISWorkspaceTypeMetadata(),
    description=_(u'Updated output workspace.'),
    arcGISDisplayName=_(u'Updated output workspace'),
    arcGISParameterDependencies=[u'outputWorkspace'])

# Public method: OceanColorLevel3SMITimeSeries.CreateClimatologicalArcGISRasters

AddMethodMetadata(OceanColorLevel3SMITimeSeries.CreateClimatologicalArcGISRasters,
    shortDescription=_(u'Creates climatological rasters for a Level 3 Standard Mapped Image (SMI) product published by the NASA GSFC OceanColor Group.'),
    longDescription=_(
u"""This tool produces rasters showing the climatological average
value (or other statistic) of a time series of Level 3 SMI images.
Given a specification of the desired Level 3 SMI time series, a
statistic, and a climatological bin definition, this tool downloads
the images, classifies them into bins, and produces a single raster
for each bin. Each cell of the raster is produced by calculating the
statistic on the values of that cell extracted from all of the rasters
in the bin.

""") + _OceanColorLevel3SMI_LongDescription % {u'name': 'tool'},
    isExposedToPythonCallers=True,
    isExposedByCOM=True,
    isExposedAsArcGISTool=True,
    arcGISDisplayName=_(u'Create Climatological Rasters for NASA OceanColor L3 SMI Product'),
    arcGISToolCategory=_(u'Data Products\\NASA GSFC OceanColor Group\\L3 Products'),
    dependencies=[ArcGISDependency(9, 1), PythonAggregatedModuleDependency('numpy')])

CopyArgumentMetadata(OceanColorLevel3SMITimeSeries.CreateArcGISRasters, u'cls', OceanColorLevel3SMITimeSeries.CreateClimatologicalArcGISRasters, u'cls')
CopyArgumentMetadata(OceanColorLevel3SMITimeSeries.CreateArcGISRasters, u'sensor', OceanColorLevel3SMITimeSeries.CreateClimatologicalArcGISRasters, u'sensor')
CopyArgumentMetadata(OceanColorLevel3SMITimeSeries.CreateArcGISRasters, u'temporalResolution', OceanColorLevel3SMITimeSeries.CreateClimatologicalArcGISRasters, u'temporalResolution')
CopyArgumentMetadata(OceanColorLevel3SMITimeSeries.CreateArcGISRasters, u'spatialResolution', OceanColorLevel3SMITimeSeries.CreateClimatologicalArcGISRasters, u'spatialResolution')
CopyArgumentMetadata(OceanColorLevel3SMITimeSeries.CreateArcGISRasters, u'product', OceanColorLevel3SMITimeSeries.CreateClimatologicalArcGISRasters, u'product')

AddArgumentMetadata(OceanColorLevel3SMITimeSeries.CreateClimatologicalArcGISRasters, u'statistic',
    typeMetadata=UnicodeStringTypeMetadata(allowedValues=[u'Count', u'Maximum', u'Mean', u'Minimum', u'Range', u'Standard Deviation', u'Sum'], makeLowercase=True),
    description=_(
u"""Statistic to calculate for each cell, one of:

* Count - number of images in which the cell had data.

* Maximum - maximum value for the cell.

* Mean - mean value for the cell, calculated as the sum divided by the
  count.

* Minimum - minimum value for the cell.

* Range - range for the cell, calculated as the maximum minus the
  minimum.

* Standard Deviation - sample standard deviation for the cell
  (i.e. the standard deviation estimated using Bessel's correction).
  In order to calculate this, there must be at least two images with
  data for the cell.

* Sum - the sum for the cell.
"""),
    arcGISDisplayName=_(u'Statistic'))

AddArgumentMetadata(OceanColorLevel3SMITimeSeries.CreateClimatologicalArcGISRasters, u'binType',
    typeMetadata=UnicodeStringTypeMetadata(allowedValues=[u'Daily', u'Monthly', u'Cumulative', u'ENSO Daily', u'ENSO Monthly', u'ENSO Cumulative'], makeLowercase=True),
    description=_(
u"""Climatology bins to use, one of:

* Daily - daily bins. Images will be classified into bins according to
  their days of the year. The number of days in each bin is determined
  by the Climatology Bin Duration parameter (which defaults to 1). The
  number of bins is calculated by dividing 365 by the bin duration. If
  there is no remainder, then that number of bins will be created;
  images for the 366th day of leap years will be counted in the bin
  that includes day 365. For example, if the bin duration is 5, 73
  bins will be created. The first will be for days 1-5, the second
  will be for days 5-10, and so on; the 73rd bin will be for days
  361-365 during normal years and 361-366 during leap years. If
  dividing 365 by the bin duration does yield a remainder, then one
  additional bin will be created to hold the remaining days. For
  example, if the bin duration is 8, 46 bins will be created. The
  first will be for days 1-8, the second for days 9-16, and so on; the
  46th will be for days 361-365 during normal years and 361-366 during
  leap years.

* Monthly - monthly bins. Images will be classified into bins according to
  their months of the year. The number of months in each bin is
  determined by the Climatology Bin Duration parameter (which defaults
  to 1). The number of bins is calculated by dividing 12 by the bin
  duration. If there is no remainder, then that number of bins will be
  created. For example, if the bin duration is 3, there will be four
  bins: January-March, April-June, July-September, and
  October-December. If there is a remainder, then one additional bin
  will be created. For example, if the bin duration is 5, 3 bins will
  be created: January-May, June-October, November-December.

* Cumulative - one bin. A single climatology raster will be calculated
  from the entire dataset. The Bin Duration parameter is ignored.

* ENSO Daily, ENSO Monthly, ENSO Cumulative - the same as above,
  except each of the bins above will be split into three, based on the
  phase of the `El Nino Southern Oscillation (ENSO) <http://en.wikipedia.org/wiki/ENSO>`_,
  as determined by the Oceanic Nino Index (ONI) calculated by the
  `NOAA NCEP Climate Prediction Center <http://www.cpc.ncep.noaa.gov/products/analysis_monitoring/ensostuff/ensoyears.shtml>`_.
  The ONI classifies each month into one of three phases: neutral, El
  Nino, or La Nina. This tool first classifies input images according
  to their dates into ENSO phases (it downloads ONI data from the
  `NOAA Earth System Research Laboratory <http://www.esrl.noaa.gov/psd/data/climateindices/list/>`_),
  then produces a climatology bin for each phase. For example, if you
  request ENSO Cumulative bins, three bins will be produced: one for
  all images occurring in neutral months, one for all in El Nino
  months, and one for all in La Nina months. If you request ENSO
  Monthly bins, 36 bins will be produced: one for each combination of
  the 12 months and the three ENSO phases.

For Daily and Monthly, to adjust when the bins start (e.g. to center a
4-bin seasonal climatology on solstices and equinoxes), use the Start
Climatology At This Day Of The Year parameter."""),
    arcGISDisplayName=_(u'Climatology bin type'))

CopyArgumentMetadata(OceanColorLevel3SMITimeSeries.CreateArcGISRasters, u'outputWorkspace', OceanColorLevel3SMITimeSeries.CreateClimatologicalArcGISRasters, u'outputWorkspace')

AddArgumentMetadata(OceanColorLevel3SMITimeSeries.CreateClimatologicalArcGISRasters, u'mode',
    typeMetadata=UnicodeStringTypeMetadata(allowedValues=[u'Add', u'Replace'], makeLowercase=True),
    description=_(
u"""Overwrite mode, one of:

* Add - create rasters that do not exist and skip those that already
  exist. This is the default.

* Replace - create rasters that do not exist and overwrite those that
  already exist. Choose this option when you want to regenerate the
  climatologies using the latest satellite images.

The ArcGIS Overwrite Outputs geoprocessing setting has no effect on
this tool. If 'Replace' is selected the rasters will be overwritten,
regardless of the ArcGIS Overwrite Outputs setting."""),
    arcGISDisplayName=_(u'Overwrite mode'))

AddArgumentMetadata(OceanColorLevel3SMITimeSeries.CreateClimatologicalArcGISRasters, u'rasterNameExpressions',
    typeMetadata=ListTypeMetadata(elementType=UnicodeStringTypeMetadata(), minLength=1),
    description=_(
u"""List of expressions specifying how the output rasters should be
named.

The default expression assumes you are storing rasters in a file
system directory and creates them in a tree structure. When storing
rasters in a directory, the final expression specifies the file name
of the raster and any preceding expressions specify subdirectories.
The extension of the final expression determines the output raster
format: .asc for ArcInfo ASCII Grid, .bmp for BMP, .gif for GIF, .img
for an ERDAS IMAGINE file, .jpg for JPEG, .jp2 for JPEG 2000, .png for
PNG, .tif for GeoTIFF, or no extension for ArcInfo Binary Grid. The
default expression uses .img.

When storing rasters in a geodatabase, you should provide only one
expression. That expression specifies the raster's name.

Each expression may contain any sequence of characters permitted by
the output workspace. Each expression may optionally contain one or
more of the following case-sensitive codes. The tool replaces the
codes with appropriate values when creating each raster:

* %(Sensor)s - name of the sensor, either "aqua", "czcs", "octs",
  "meris", "seawifs", or "terra".

* %(SensorCode)s - abbreviation for the sensor, either "A", "C", "O",
  "S", or "T", used in NASA's file naming scheme.

* %(TemporalResolution)s - temporal resolution, either "daily",
  "8day", "monthly", or "annual".

* %(TemporalResolutionCode)s - abbreviation for the temporal
  resolution, either "DAY", "8D", "MO", or "YR", used in NASA's file
  naming scheme.

* %(SpatialResolution)s - spatial resolution, either "4km" or "9km".

* %(SpatialResolutionCode)s - abbreviation for the spatial
  resolution, either "4" or "9", used in NASA's file naming scheme for
  certain products.
  
* %(ProductCode)s - product code for the Level 3 Standard Mapped Image
  (SMI) product represented in the output raster.

* %(ClimatologyBinType)s - type of the climatology bin, either "Daily"
  if 1-day bins, "Xday" if multi-day bins (X is replaced by the
  duration), "Monthly" if 1-month bins, "Xmonth" if multi-month bins,
  or "Cumulative". If an ENSO bin type is used, "ENSO\\_" will be
  prepended to those strings (e.g. "ENSO_Daily", "ENSO_Monthly").

* %(ClimatologyBinName)s - name of the climatology bin corresponding
  represented by the output raster, either "dayXXX" for 1-day bins
  (XXX is replaced by the day of the year), "daysXXXtoYYY" for
  multi-day bins (XXX is replaced by the first day of the bin, YYY is
  replaced by the last day), "monthXX" for 1-month bins (XX is
  replaced by the month), "monthXXtoYY" (XX is replaced by the first
  month of the bin, YY by the last month), or "cumulative". If an ENSO
  bin type is used, "neutral\\_", "ElNino\\_", and "LaNina\\_" will be
  prepended to those strings for each of the three ENSO phased rasters
  (e.g. "neutral_cumulative", "ElNino_cumulative", and
  "LaNina_cumulative" when "ENSO Cumulative" bins are requested).

* %(Statistic)s - statistic that was calculated, in lowercase and with
  spaces replaced by underscores; one of: "count", "maximum", "mean",
  "minimum", "range", "standard_deviation", "Sum".

If the Bin Type is "Daily", the following additional codes are
available:

* %(FirstDay)i - first day of the year of the climatology bin
  represented by the output raster.

* %(LastDay)i - last day of the year of the climatology bin
  represented by the output raster. For 1-day climatologies, this will
  be the same as %(FirstDay)i.

If the Bin Type is "Monthly", the following additional codes are
available:

* %(FirstMonth)i - first month of the climatology bin represented by
  the output raster.

* %(DayOfFirstMonth)i - first day of the first month of the
  climatology bin represented by the output raster.

* %(LastMonth)i - last month of the climatology bin represented by
  the output raster.

* %(DayOfLastMonth)i - last day of the last month of the climatology
  bin represented by the output raster.

Note that the additional codes are integers and may be formatted using
"printf"-style formatting codes. For example, to format the FirstDay
as a three-digit number with leading zeros::

    %(FirstDay)03i
"""),
    arcGISDisplayName=_(u'Raster name expressions'))

CopyArgumentMetadata(OceanColorLevel3SMITimeSeries.CreateArcGISRasters, u'cacheDirectory', OceanColorLevel3SMITimeSeries.CreateClimatologicalArcGISRasters, u'cacheDirectory')

AddArgumentMetadata(OceanColorLevel3SMITimeSeries.CreateClimatologicalArcGISRasters, u'binDuration',
    typeMetadata=IntegerTypeMetadata(minValue=1),
    description=_(
u"""Duration of each bin, in days or months, when the Bin Type is
Daily or Monthly, respectively. The default is 1. See the Bin Type
parameter for more information."""),
    arcGISDisplayName=_(u'Climatology bin duration'),
    arcGISCategory=_(u'Climatology options'))

AddArgumentMetadata(OceanColorLevel3SMITimeSeries.CreateClimatologicalArcGISRasters, u'startDayOfYear',
    typeMetadata=IntegerTypeMetadata(minValue=1, maxValue=365),
    description=_(
u"""Use this parameter to create bin defintions that deviate from the
traditional calendar. The interpretation of this parameter depends on
the Bin Type:

* Daily - this parameter defines the day of the year of the first
  climatology bin. For example, if this parameter is 100 and the Bin
  Duration is 10, the first bin will be numbered 100-109. The bin
  spanning the end of the year will be numbered 360-004. The last bin
  will be numbered 095-099. To define a four-bin climatology with bins
  that are centered approximately on the equinoxes and solstices
  (i.e., a seasonal climatology), set the Bin Duration to 91 and the
  start day to 36 (February 5). This will produce bins with dates
  036-126, 127-217, 218-308, and 309-035.

* Monthly - this parameter defines the day of the year of the first
  climatology bin, and the day of the month of that bin will be used
  as the first day of the month of all of the bins. For example, if
  this parameter is 46, which is February 15, and the Bin Duration is
  1, then the bins will be February 15 - March 14, March 15 - April
  14, April 15 - May 14, and so on. Calculations involving this
  parameter always assume a 365 day year (a non-leap year). To define
  a four-bin climatology using the months traditionally associated
  with spring, summer, fall, and winter in many northern hemisphere
  cultures, set the Bin Duration to 3 and the start day to 60 (March
  1). This will produce bins with months 03-05, 06-08, 09-11, and
  12-02.

* Cumulative - this parameter is ignored.
"""),
    arcGISDisplayName=_(u'Start climatology at this day of the year'),
    arcGISCategory=_(u'Climatology options'))

CopyArgumentMetadata(OceanColorLevel3SMITimeSeries.CreateArcGISRasters, u'rotationOffset', OceanColorLevel3SMITimeSeries.CreateClimatologicalArcGISRasters, u'rotationOffset')
CopyArgumentMetadata(OceanColorLevel3SMITimeSeries.CreateArcGISRasters, u'spatialExtent', OceanColorLevel3SMITimeSeries.CreateClimatologicalArcGISRasters, u'spatialExtent')
CopyArgumentMetadata(OceanColorLevel3SMITimeSeries.CreateArcGISRasters, u'startDate', OceanColorLevel3SMITimeSeries.CreateClimatologicalArcGISRasters, u'startDate')
CopyArgumentMetadata(OceanColorLevel3SMITimeSeries.CreateArcGISRasters, u'endDate', OceanColorLevel3SMITimeSeries.CreateClimatologicalArcGISRasters, u'endDate')
CopyArgumentMetadata(OceanColorLevel3SMITimeSeries.CreateArcGISRasters, u'timeout', OceanColorLevel3SMITimeSeries.CreateClimatologicalArcGISRasters, u'timeout')
CopyArgumentMetadata(OceanColorLevel3SMITimeSeries.CreateArcGISRasters, u'maxRetryTime', OceanColorLevel3SMITimeSeries.CreateClimatologicalArcGISRasters, u'maxRetryTime')
CopyArgumentMetadata(OceanColorLevel3SMITimeSeries.CreateArcGISRasters, u'calculateStatistics', OceanColorLevel3SMITimeSeries.CreateClimatologicalArcGISRasters, u'calculateStatistics')
CopyArgumentMetadata(OceanColorLevel3SMITimeSeries.CreateArcGISRasters, u'buildPyramids', OceanColorLevel3SMITimeSeries.CreateClimatologicalArcGISRasters, u'buildPyramids')

CopyResultMetadata(OceanColorLevel3SMITimeSeries.CreateArcGISRasters, u'updatedOutputWorkspace', OceanColorLevel3SMITimeSeries.CreateClimatologicalArcGISRasters, u'updatedOutputWorkspace')

# Public method: OceanColorLevel3SMITimeSeries.InterpolateAtArcGISPoints

AddMethodMetadata(OceanColorLevel3SMITimeSeries.InterpolateAtArcGISPoints,
    shortDescription=_(u'Interpolates the values of a Level 3 Standard Mapped Image (SMI) product published by the NASA GSFC OceanColor Group at points.'),
    longDescription=_(
u"""Given a sensor name, temporal resolution, spatial resolution, and
desired Level 3 SMI product, this tool interpolates the value of that
product at the given points. This tool performs the same basic
operation as the ArcGIS Spatial Analyst's Extract Values to Points
tool, but it downloads and reads HDF files from NASA's servers rather
than reading rasters stored on your machine.

""") + _OceanColorLevel3SMI_LongDescription % {u'name': 'tool'},
    isExposedToPythonCallers=True,
    isExposedByCOM=True,
    isExposedAsArcGISTool=True,
    arcGISDisplayName=_(u'Interpolate NASA OceanColor L3 SMI Product at Points'),
    arcGISToolCategory=_(u'Data Products\\NASA GSFC OceanColor Group\\L3 Products'),
    dependencies=[ArcGISDependency(9, 1, requiresCOMInstantiation=True), PythonAggregatedModuleDependency('numpy')])

CopyArgumentMetadata(OceanColorLevel3SMITimeSeries.CreateArcGISRasters, u'cls', OceanColorLevel3SMITimeSeries.InterpolateAtArcGISPoints, u'cls')
CopyArgumentMetadata(OceanColorLevel3SMITimeSeries.__init__, u'sensor', OceanColorLevel3SMITimeSeries.InterpolateAtArcGISPoints, u'sensor')
CopyArgumentMetadata(OceanColorLevel3SMITimeSeries.__init__, u'temporalResolution', OceanColorLevel3SMITimeSeries.InterpolateAtArcGISPoints, u'temporalResolution')
CopyArgumentMetadata(OceanColorLevel3SMITimeSeries.__init__, u'spatialResolution', OceanColorLevel3SMITimeSeries.InterpolateAtArcGISPoints, u'spatialResolution')
CopyArgumentMetadata(OceanColorLevel3SMITimeSeries.__init__, u'product', OceanColorLevel3SMITimeSeries.InterpolateAtArcGISPoints, u'product')

AddArgumentMetadata(OceanColorLevel3SMITimeSeries.InterpolateAtArcGISPoints, u'points',
    typeMetadata=ArcGISFeatureLayerTypeMetadata(mustExist=True, allowedShapeTypes=[u'Point']),
    description=_(
u"""Points at which values should be interpolated.

OceanColor images use a geographic coordinate system with the WGS 1984
datum. It is recommended but not required that the points use the same
coordinate system. If they do not, this tool will attempt to project
the points to the OceanColor coordinate system prior to doing the
interpolation. This may fail if a datum transformation is required, in
which case you will have to manually project the points to the
OceanColor coordinate system before using this tool."""),
    arcGISDisplayName=_(u'Point features'))

AddArgumentMetadata(OceanColorLevel3SMITimeSeries.InterpolateAtArcGISPoints, u'valueField',
    typeMetadata=ArcGISFieldTypeMetadata(mustExist=True, allowedFieldTypes=[u'short', u'long', u'float', u'double']),
    description=_(
u"""Field of the points to receive the interpolated values.

The field must have a floating-point or integer data type. If the
field cannot represent the interpolated value at full precision, the
closest approximation will be stored and a warning will be issued.
This will happen, for example, when you interpolate values into an
integer field."""),
    arcGISDisplayName=_(u'Field to receive the interpolated values'),
    arcGISParameterDependencies=[u'points'])

AddArgumentMetadata(OceanColorLevel3SMITimeSeries.InterpolateAtArcGISPoints, u'tField',
    typeMetadata=ArcGISFieldTypeMetadata(mustExist=True, allowedFieldTypes=[u'date']),
    description=_(
u"""Field of the points that specifies the date and time of the point.

The field must have a date or datetime data type. If the field can
only represent dates with no time component, the time will assumed to
be 00:00:00."""),
    arcGISDisplayName=_(u'Date field'),
    arcGISParameterDependencies=[u'points'])

AddArgumentMetadata(OceanColorLevel3SMITimeSeries.InterpolateAtArcGISPoints, u'method',
    typeMetadata=UnicodeStringTypeMetadata(allowedValues=[u'Nearest', u'Linear'], makeLowercase=True),
    description=_(
u"""Interpolation method to use, one of:

* Nearest - nearest neighbor interpolation. The interpolated value
  will simply be the value of the cell that contains the point. This
  is the default.

* Linear - linear interpolation (also known as trilinear
  interpolation). This method averages the values of the eight nearest
  cells in the x, y, and time dimensions, weighting the contribution
  of each cell by the area of it that would be covered by a
  hypothetical cell centered on the point being interpolated. If the
  cell containing the point contains NoData, the result is NoData. If
  any of the other seven cells contain NoData, they are omitted from
  the average, and the result is based on the weighted average of the
  cells that do contain data. This is the same algorithm implemented
  by the ArcGIS Spatial Analyst's Extract Values to Points tool.
"""),
    arcGISDisplayName=_(u'Interpolation method'))

CopyArgumentMetadata(OceanColorLevel3SMITimeSeries.__init__, u'cacheDirectory', OceanColorLevel3SMITimeSeries.InterpolateAtArcGISPoints, u'cacheDirectory')

AddArgumentMetadata(OceanColorLevel3SMITimeSeries.InterpolateAtArcGISPoints, u'where',
    typeMetadata=SQLWhereClauseTypeMetadata(canBeNone=True),
    description=_(
u"""SQL WHERE clause expression that specifies the subset of points to
use. If this parameter is not provided, all of the points will be
used.

The exact syntax of this expression depends on the type of feature
class you're using. ESRI recommends you reference fields using the
following syntax:

* For shapefiles, ArcInfo coverages, or feature classes stored in file
  geodatabases, ArcSDE geodatabases, or ArcIMS, enclose field names in
  double quotes: "MY_FIELD"

* For feature classes stored in personal geodatabases, enclose field
  names in square brackets: [MY_FIELD].
"""),
    arcGISDisplayName=_(u'Where clause'),
    arcGISCategory=_(u'Interpolation options'),
    arcGISParameterDependencies=[u'points'])

AddArgumentMetadata(OceanColorLevel3SMITimeSeries.InterpolateAtArcGISPoints, u'noDataValue',
    typeMetadata=FloatTypeMetadata(canBeNone=True),
    description=_(
u"""Value to use when the interpolated value is NoData.

If a value is not provided for this parameter, a database NULL value
will be stored in the field when the interpolated value is NoData. If
the field cannot store NULL values, as is the case with shapefiles,
the value -9999 will be used."""),
    arcGISDisplayName=_(u'Value to use when the interpolated value is NoData'),
    arcGISCategory=_(u'Interpolation options'))

CopyArgumentMetadata(OceanColorLevel3SMITimeSeries.__init__, u'timeout', OceanColorLevel3SMITimeSeries.InterpolateAtArcGISPoints, u'timeout')
CopyArgumentMetadata(OceanColorLevel3SMITimeSeries.__init__, u'maxRetryTime', OceanColorLevel3SMITimeSeries.InterpolateAtArcGISPoints, u'maxRetryTime')

AddArgumentMetadata(OceanColorLevel3SMITimeSeries.InterpolateAtArcGISPoints, u'orderByFields',
    typeMetadata=ListTypeMetadata(elementType=ArcGISFieldTypeMetadata(mustExist=True), minLength=1, canBeNone=True),
    description=_(
u"""Fields for defining the order in which the points are processed.

The points may be processed faster if they are ordered
spatiotemporally, such that points that are close in space and time
are processed sequentially. Ordering the points this way increases the
probability that the value of a given point can be interpolated from
data that is cached in memory, rather than from data that must be read
from the disk or network, which is much slower. Choose fields that
faciliate this. For example, if your points represent the locations of
animals tracked by satellite telemetry, order the processing first by
the animal ID and then by the transmission date or number.

If you omit this parameter, the Date Field will be used automatically.

This parameter requires ArcGIS 9.2 or later."""),
    arcGISDisplayName=_(u'Order by fields'),
    arcGISCategory=_(u'Performance tuning options'),
    arcGISParameterDependencies=[u'points'],
    dependencies=[ArcGISDependency(9, 2)])

AddArgumentMetadata(OceanColorLevel3SMITimeSeries.InterpolateAtArcGISPoints, u'numBlocksToCacheInMemory',
    typeMetadata=IntegerTypeMetadata(minValue=0, canBeNone=True),
    description=_(
u"""Maximum number of blocks of data to cache in memory.

To minimize the number of times that the disk or network must be
accessed, this tool employs a simple caching strategy, in addition to
disk caching described by the Cache Directory parameter. When it
processes the first point, it reads a square block of cells centered
on that point and caches it in memory. When it processes the second
and subsequent points, it first checks whether the cells needed for
that point are contained by the block cached in memory. If so, it
processes that point using the in-memory block, rather than reading
from disk or the network again. If not, it reads another square block
centered on that point and adds it to the cache.

The tool processes the remaining points, adding additional blocks to
the cache, as needed. To prevent the cache from exhausing all memory,
it is only permitted to grow to the size specified by this parameter.
When the cache is full but a new block is needed, the oldest block is
discarded to make room for the newest block.

The maximum size of the cache in bytes may be calculated by
multiplying this parameter by 4 and by the block size parameters. For
example, if this parameter is 128 and the blocks are x=32 by y=32 by
t=2, the maximum size of the cache is 1048576 bytes (1 MB).

If this parameter is 0, no blocks will be cached in memory."""),
    arcGISDisplayName=_(u'Number of blocks of data to cache in memory'),
    arcGISCategory=_(u'Performance tuning options'))

AddArgumentMetadata(OceanColorLevel3SMITimeSeries.InterpolateAtArcGISPoints, u'xBlockSize',
    typeMetadata=IntegerTypeMetadata(minValue=0, canBeNone=True),
    description=_(
u"""Size of the blocks of data to cache in memory, in the x direction
(longitude). The size is given as the number of cells.

If this parameter is 0, no blocks will be cached in memory."""),
    arcGISDisplayName=_(u'In-memory cache block size, in X direction'),
    arcGISCategory=_(u'Performance tuning options'))

AddArgumentMetadata(OceanColorLevel3SMITimeSeries.InterpolateAtArcGISPoints, u'yBlockSize',
    typeMetadata=IntegerTypeMetadata(minValue=0, canBeNone=True),
    description=_(
u"""Size of the blocks of data to cache in memory, in the y direction
(latitude). The size is given as the number of cells.

If this parameter is 0, no blocks will be cached in memory."""),
    arcGISDisplayName=_(u'In-memory cache block size, in Y direction'),
    arcGISCategory=_(u'Performance tuning options'))

AddArgumentMetadata(OceanColorLevel3SMITimeSeries.InterpolateAtArcGISPoints, u'tBlockSize',
    typeMetadata=IntegerTypeMetadata(minValue=0, canBeNone=True),
    description=_(
u"""Size of the blocks of data to cache in memory, in the t direction
(time). The size is given as the number of cells.

If this parameter is 0, no blocks will be cached in memory."""),
    arcGISDisplayName=_(u'In-memory cache block size, in T direction'),
    arcGISCategory=_(u'Performance tuning options'))

AddResultMetadata(OceanColorLevel3SMITimeSeries.InterpolateAtArcGISPoints, u'updatedPoints',
    typeMetadata=ArcGISFeatureLayerTypeMetadata(),
    description=_(u'Updated points.'),
    arcGISDisplayName=_(u'Updated points'),
    arcGISParameterDependencies=[u'points'])

###############################################################################
# Metadata: OceanColorLevel2Converter class
###############################################################################

AddClassMetadata(OceanColorLevel2Converter,
    shortDescription=_(u'This class contains methods for converting L2 files (swath granules) published by NASA OceanColor to ArcGIS-compatible formats.'))

# Public method: OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints

AddMethodMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints,
    shortDescription=_(u'Converts a NASA OceanColor MODIS L2 LAC SST file to an ArcGIS point feature class.'),
    longDescription=_(u'TODO: Write long description.'),
    isExposedToPythonCallers=True,
    isExposedByCOM=True,
    isExposedAsArcGISTool=True,
    arcGISDisplayName=_(u'Convert MODIS L2 LAC SST to ArcGIS Points'),
    arcGISToolCategory=_(u'Data Products\\NASA GSFC OceanColor Group\\L2 Products'),
    dependencies=[ArcGISDependency(9, 1), PythonAggregatedModuleDependency('numpy')])

AddArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'cls',
    typeMetadata=ClassOrClassInstanceTypeMetadata(cls=OceanColorLevel2Converter),
    description=_(u'OceanColorLevel2Converter class or instance.'))

AddArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'hdf',
    typeMetadata=FileTypeMetadata(mayBeCompressed=True, mustExist=True),
    description=_(
u"""MODIS L2 LAC SST file to convert.

MODIS L2 LAC SST files have the extension .L2_LAC_SST.bz2 or
.L2_LAC_SST4.bz2 and can be downloaded from the
`OceanColor L2 Browser <http://oceancolor.gsfc.nasa.gov/cgi/browse.pl>`_.
You do not need to decompress the files; they will be decompressed
automatically.

This tool only works on the MODIS L2 LAC SST product; it will not work
on other L2 products."""),
    arcGISDisplayName=_(u'L2 LAC SST file'))

AddArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'points',
    typeMetadata=ArcGISFeatureClassTypeMetadata(mustBeDifferentThanArguments=[u'hdf'], deleteIfParameterIsTrue=u'overwriteExisting', createParentDirectories=True),
    description=_(
u"""Output point feature class to create.

Each pixel that is not masked will be converted to a point. The points
will have an attribute named sst, qual_sst, or l2_flags that contains
the pixel values. The next parameter of this tool determines which of
these attributes will be used."""),
    direction = u'Output',
    arcGISDisplayName=_(u'Output point feature class'))

AddArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'variable',
    typeMetadata=UnicodeStringTypeMetadata(allowedValues=[u'sst', u'qual_sst', u'l2_flags'], makeLowercase=True),
    description=_(
u"""Variable to extract from the file, one of:

* sst - sea surface temperature, in degrees Celsius.

* qual_sst - quality level of the SST estimate, with values 0 (best),
  1 (good), 2 (questionable), 3 (bad), or 4 (failed or not computed).
  `This page <http://oceancolor.gsfc.nasa.gov/DOCS/modis_sst/>`_ gives
  the details of how the quality levels are computed.

* l2_flags - level-2 flags for the pixel. These are the results of various
  quality tests that are mostly intended to be used with chlorophyll
  and not SST, but because they are always present in the SST file, we
  give you the option of extracting them. The page
  `here <http://oceancolor.gsfc.nasa.gov/VALIDATION/flags.html>`_
  documents the level-2 flags.
"""),
    arcGISDisplayName=_(u'Variable'))

AddArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'spatialExtent',
    typeMetadata=EnvelopeTypeMetadata(canBeNone=True),
    description=_(
u"""Geographic bounding box of pixels to include. Pixels falling
outside of this bounding box will be ignored and excluded from the
output. If no bounding box is provided, all pixels in the file will be
used."""),
    arcGISDisplayName=_(u'Geographic extent'))

AddArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'qualityLevel',
    typeMetadata=IntegerTypeMetadata(allowedValues=[0, 1, 2, 3, 4]),
    description=_(
u"""Minimum pixel quality level to use, one of:

* 0 - best
* 1 - good (the default)
* 2 - questionable
* 3 - bad
* 4 - failed or SST not computed

All pixels with a quality level that is worse than the minimum will be
masked and excluded from the output."""),
    arcGISDisplayName=_(u'Minimum quality level'))

AddArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'recoverBadPixels',
    typeMetadata=BooleanTypeMetadata(),
    description=_(
u"""This is an experimental option, not known to be endorsed by NASA,
designed to make MODIS L2 SST more usable for SST front detection. If
this option is enabled, it will only be used for daytime files for the
Aqua satellite. It is ignored for nighttime files and for the Terra
satellite.

If this option is enabled, the tool will use the
`l2_flags <http://oceancolor.gsfc.nasa.gov/VALIDATION/flags.html>`_
computed for the chlorophyll product to try to identify pixels flagged
with SST quality levels 1, 2, or 3 that should have been assigned a
quality level of 0 but were not, presumably because they occurred
along a frontal boundary.

For each pixel flagged with quality level 1, 2, or 3, the tool
examines the l2_flags. The parameters following this one specify which
l2_flags are examined. If none of those flags are set, the tool treats
the pixel as if it had a quality level of 0, effectively "unmasking"
it. The philosophy here is: if the quality tests designed for
chlorophyll did not indicate land, cloud, ice or some other major
problem, assume the SST pixel was mistakenly marked with a low quality
level due to it being along a front. NASA acknowledges `here
<http://oceancolor.gsfc.nasa.gov/DOCS/modis_sst/>`_ (under the heading
Quality Tests) that frontal boundaries may cause SST quality tests
such as SSTREFDIFF, SSTREFVDIFF, BTNONUNIF, and BTVNONUNIF to fail,
leading to a quality level of 2 or 3 being assigned.

This "unmasking" technique cannot be applied to nighttime files
because the required l2_flags are only computed properly for daytime
files. The color of the ocean can only be observed when the sun
illuminates the ocean. We assumed that it also cannot be succesfully
applied to Terra files because the MODIS instrument on Terra has
defects that prevent it from estimating chlorophyll properly. For more
information about those defects, see
`this message <http://oceancolor.gsfc.nasa.gov/forum/oceancolor/topic_show.pl?tid=3734>`_.

We have configured the remaining parameters of this tool to use the
l2_flags that we believe provide the best compromise between
maximizing the number of pixels that are unmasked and minimizing the
number that should have remained masked. If you do not like the
results, we encourage you to experiment with different combinations of
l2_flags. Under the documentation for each parameter, we provide a
justification for our default choice for that L2 flag."""),
    arcGISDisplayName=_(u'Use l2_flags to unmask frontal boundaries in Aqua daytime images'))

AddArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkATMFAIL',
    typeMetadata=BooleanTypeMetadata(),
    description=_(
u"""If both this option and the option to use l2_flags to unmask
daytime images are enabled, then the ATMFAIL flag will be checked, and
if that flag is set, the pixel will remain masked.

This option is enabled by default.

`Here <http://oceancolor.gsfc.nasa.gov/VALIDATION/flags.html>`_, NASA
describes the ATMFAIL flag as "atmospheric correction failure". We did
not research the circumstances under which this failure might occur;
the details might appear in
'this document <http://oceancolor.gsfc.nasa.gov/SeaWiFS/TECH_REPORTS/PLVol22.pdf>_.
We assumed that if the atmospheric correction algorithm failed, it was
likely that the atmosphere was in an anomalous or unexpected state, or
that sufficient data were not available to correct for the state of
atmosphere, and that the SST estimates would also have significant
error. This may be a bad assumption if the atmospheric correction
algorithm is designed explicitly for chlorophyll but we elected to be
conservative and leave the pixels masked."""),
    arcGISDisplayName=_(u'ATMFAIL - Atmospheric correction failure'),
    arcGISCategory=_(u'l2_flags to check'))

AddArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkLAND',
    typeMetadata=BooleanTypeMetadata(),
    description=_(
u"""If both this option and the option to use l2_flags to unmask
daytime images are enabled, then the LAND flag will be checked, and if
that flag is set, the pixel will remain masked.

This option is enabled by default.

`Here <http://oceancolor.gsfc.nasa.gov/VALIDATION/flags.html>`_, NASA
describes the LAND flag as "pixel is over land".
"""),
    arcGISDisplayName=_(u'LAND - Pixel is over land'),
    arcGISCategory=_(u'l2_flags to check'))

AddArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkPRODWARN',
    typeMetadata=BooleanTypeMetadata(),
    description=_(
u"""If both this option and the option to use l2_flags to unmask
daytime images are enabled, then the PRODWARN flag will be checked,
and if that flag is set, the pixel will remain masked.

This option is disabled by default.

`Here <http://oceancolor.gsfc.nasa.gov/VALIDATION/flags.html>`_, NASA
describes the PRODWARN flag as "one or more product warnings". We do
not know what this means, but assume it is derived from the other
flags, or is specific to the ocean color algorithms. In either case,
we believed it was not appropriate to check this flag for SST
pixels."""),
    arcGISDisplayName=_(u'PRODWARN - One or more product warnings'),
    arcGISCategory=_(u'l2_flags to check'))

AddArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkHIGLINT',
    typeMetadata=BooleanTypeMetadata(),
    description=_(
u"""If both this option and the option to use l2_flags to unmask
daytime images are enabled, then the HIGLINT flag will be checked, and
if that flag is set, the pixel will remain masked.

This option is disabled by default.

`Here <http://oceancolor.gsfc.nasa.gov/VALIDATION/flags.html>`_, NASA
describes the HIGLINT flag as "high sun glint" and says "this bit is
set if glint reflectance exceeds 0.005". As we understand it, the
longwave SST algorithm used for daytime SST images is not affected by
sun glint, so we disabled the checking of this flag by default.
For an image that illustrates this, see page 14 of
`this presentation <http://oceancolor.gsfc.nasa.gov/DOCS/Presentations/seadas_umbc_franz.pdf>`_.
"""),
    arcGISDisplayName=_(u'HIGLINT - High sun glint'),
    arcGISCategory=_(u'l2_flags to check'))

AddArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkHILT',
    typeMetadata=BooleanTypeMetadata(),
    description=_(
u"""If both this option and the option to use l2_flags to unmask
daytime images are enabled, then the HILT flag will be checked, and if
that flag is set, the pixel will remain masked.

This option is enabled by default.

`Here <http://oceancolor.gsfc.nasa.gov/VALIDATION/flags.html>`_, NASA
describes the HILT flag as "observed radiance very high or saturated"
and says "this bit is set if digital count value in bands 7 and 8 is
above knee".

TODO: Test this option to determine an appropriate default. """),
    arcGISDisplayName=_(u'HILT - Observed radiance very high or saturated'),
    arcGISCategory=_(u'l2_flags to check'))

AddArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkHISATZEN',
    typeMetadata=BooleanTypeMetadata(),
    description=_(
u"""If both this option and the option to use l2_flags to unmask
daytime images are enabled, then the HISATZEN flag will be checked, and
if that flag is set, the pixel will remain masked.

This option is enabled by default.

`Here <http://oceancolor.gsfc.nasa.gov/VALIDATION/flags.html>`_, NASA
describes the HISATZEN flag as "high sensor view zenith angle" and
says "this bit is set if 60 degree satellite zenith threshold is
reached". NASA also checks the sensor zenith angle during 
`SST processing <http://oceancolor.gsfc.nasa.gov/DOCS/modis_sst/>`_.
If the angle exceeds 55 degrees, the SST quality is set to 1 (good);
if it exceeds 75 degrees, it is set to 3 (bad).

We defer to NASA's finding that an angle of 75 degrees or greater
results in bad SST estimates. Therefore we enabled checking the
HISATZEN chlorophyll flag by default. As a result, when the SST
quality level is 2, 3, or 4, and the zenith angle is greater than 60,
the pixel will not be unmasked. This is not a perfect solution, but it
errs on conservative side."""),
    arcGISDisplayName=_(u'HISATZEN - High sensor view zenith angle'),
    arcGISCategory=_(u'l2_flags to check'))

AddArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkCOASTZ',
    typeMetadata=BooleanTypeMetadata(),
    description=_(
u"""If both this option and the option to use l2_flags to unmask
daytime images are enabled, then the COASTZ flag will be checked, and
if that flag is set, the pixel will remain masked.

This option is disabled by default.

`Here <http://oceancolor.gsfc.nasa.gov/VALIDATION/flags.html>`_, NASA
describes the COASTZ flag as "pixel is in shallow water". As far as we
know, MODIS SST estimates are not problematic in shallow water,
therefore we disabled the checking of this flag by default."""),
    arcGISDisplayName=_(u'COASTZ - Pixel is in shallow water'),
    arcGISCategory=_(u'l2_flags to check'))

AddArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkSTRAYLIGHT',
    typeMetadata=BooleanTypeMetadata(),
    description=_(
u"""If both this option and the option to use l2_flags to unmask
daytime images are enabled, then the STRAYLIGHT flag will be checked, and
if that flag is set, the pixel will remain masked.

This option is enabled by default.

`Here <http://oceancolor.gsfc.nasa.gov/VALIDATION/flags.html>`_, NASA
describes the STRAYLIGHT flag as "straylight contamination is likely"
and then links to
'this document <http://oceancolor.gsfc.nasa.gov/REPROCESSING/Aqua/R1/modisa_repro1_stlight.html>`_.
In short, stray light refers to the problem in which pixels that are a
few kilometers from bright sources such as clouds or land are
erroneously bright, due to imperfect sensor optics. We decided to
enable checking the STRAYLIGHT flag by default after noticing that the
CLDICE flag did not completely mask certain clouds. Sometimes there
was a halo of cold pixels around the clouds; during SST front
detection, a front would be erroneously detected along these pixels.
When we enabled checking the STRAYLIGHT flag the clouds were
effectively buffered, masking the cold halo.

Unfortunately, this strategy is not perfect. The STRAYLIGHT flag masks
a lot of pixels that are not obviously invalid, particularly close to
shorelines. If you are working in a coastal area and require SST
estimates very close to shore, you may have to disable the STRAYLIGHT
flag to obtain any estimates."""),
    arcGISDisplayName=_(u'STRAYLIGHT - Straylight contamination is likely'),
    arcGISCategory=_(u'l2_flags to check'))

AddArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkCLDICE',
    typeMetadata=BooleanTypeMetadata(),
    description=_(
u"""If both this option and the option to use l2_flags to unmask
daytime images are enabled, then the CLDICE flag will be checked, and
if that flag is set, the pixel will remain masked.

This option is enabled by default.

`Here <http://oceancolor.gsfc.nasa.gov/VALIDATION/flags.html>`_, NASA
describes the CLDICE flag as "probable cloud or ice contamination". We
enabled checking this flag by default, under the rationale that if the
pixel appears to be contaminated by clouds or ice during chlorophyll
processing, the contamination will also affect the SST estimate."""),
    arcGISDisplayName=_(u'CLDICE - Probable cloud or ice contamination'),
    arcGISCategory=_(u'l2_flags to check'))

AddArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkCOCCOLITH',
    typeMetadata=BooleanTypeMetadata(),
    description=_(
u"""If both this option and the option to use l2_flags to unmask
daytime images are enabled, then the COCCOLITH flag will be checked, and
if that flag is set, the pixel will remain masked.

This option is disabled by default.

`Here <http://oceancolor.gsfc.nasa.gov/VALIDATION/flags.html>`_, NASA
describes the COCCOLITH flag as "coccolithophores detected". As far as
we know, MODIS SST estimates are not problematic when coccolithophores
are present, therefore we disabled the checking of this flag by
default."""),
    arcGISDisplayName=_(u'COCCOLITH - Coccolithophores detected'),
    arcGISCategory=_(u'l2_flags to check'))

AddArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkTURBIDW',
    typeMetadata=BooleanTypeMetadata(),
    description=_(
u"""If both this option and the option to use l2_flags to unmask
daytime images are enabled, then the TURBIDW flag will be checked, and
if that flag is set, the pixel will remain masked.

This option is disabled by default.

`Here <http://oceancolor.gsfc.nasa.gov/VALIDATION/flags.html>`_, NASA
describes the TURBIDW flag as "turbid water detected". As far as we
know, MODIS SST estimates are not problematic in turbid water,
therefore we disabled the checking of this flag by default."""),
    arcGISDisplayName=_(u'TURBIDW - Turbid water detected'),
    arcGISCategory=_(u'l2_flags to check'))

AddArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkHISOLZEN',
    typeMetadata=BooleanTypeMetadata(),
    description=_(
u"""If both this option and the option to use l2_flags to unmask
daytime images are enabled, then the HISOLZEN flag will be checked,
and if that flag is set, the pixel will remain masked.

This option is enabled by default.

`Here <http://oceancolor.gsfc.nasa.gov/VALIDATION/flags.html>`_, NASA
describes the HISOLZEN flag as "high solar zenith". While the MODIS
`SST processing <http://oceancolor.gsfc.nasa.gov/DOCS/modis_sst/>`_
does not appear to be negatively impacted by high solar zenith angles,
we are not certain that the relevant chlorophyll l2_flags are correct
at high solar zenith angles. Therefore, to be conservative, we enabled
checking this flag by default."""),
    arcGISDisplayName=_(u'HISOLZEN - High solar zenith'),
    arcGISCategory=_(u'l2_flags to check'))

AddArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkLOWLW',
    typeMetadata=BooleanTypeMetadata(),
    description=_(
u"""If both this option and the option to use l2_flags to unmask
daytime images are enabled, then the LOWLW flag will be checked,
and if that flag is set, the pixel will remain masked.

This option is enabled by default.

`Here <http://oceancolor.gsfc.nasa.gov/VALIDATION/flags.html>`_, NASA
describes the LOWLW flag as "very low water-leaving radiance (cloud
shadow)". We could find little NASA documentation for this flag.
Although we know of no specific reason why SST estimates would be bad
in pixels that are in the shadow of a cloud, we want to err on the
side of caution in unmasking pixels that are close to clouds, to
minimize the chances of unmasking cloudy pixels. Therefore we enabled
this option by default."""),
    arcGISDisplayName=_(u'LOWLW - Very low water-leaving radiance (cloud shadow)'),
    arcGISCategory=_(u'l2_flags to check'))

AddArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkCHLFAIL',
    typeMetadata=BooleanTypeMetadata(),
    description=_(
u"""If both this option and the option to use l2_flags to unmask
daytime images are enabled, then the CHLFAIL flag will be checked, and
if that flag is set, the pixel will remain masked.

This option is disabled by default.

`Here <http://oceancolor.gsfc.nasa.gov/VALIDATION/flags.html>`_, NASA
describes the CHLFAIL flag as "derived product algorithm failure".
Presumably the derived product is chlorophyll.
`This document <http://oceancolor.gsfc.nasa.gov/SeaWiFS/TECH_REPORTS/PLVol22.pdf>`_
describes some scenarios that will set the CHLFAIL flag. Those
scenarios seem to be very specific to chlorophyll estimation. We see
no reason why the SST estimation would be affected, therefore we 
disabled the checking of this flag by default."""),
    arcGISDisplayName=_(u'CHLFAIL - Derived product algorithm failure'),
    arcGISCategory=_(u'l2_flags to check'))

AddArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkNAVWARN',
    typeMetadata=BooleanTypeMetadata(),
    description=_(
u"""If both this option and the option to use l2_flags to unmask
daytime images are enabled, then the NAVWARN flag will be checked, and
if that flag is set, the pixel will remain masked.

This option is enabled by default.

`Here <http://oceancolor.gsfc.nasa.gov/VALIDATION/flags.html>`_, NASA
describes the NAVWARN flag as "navigation quality is reduced". In
satellite remote sensing, navigation refers to the problem of
determining accurate and precise geographic coordinates for each
measurement taken by the sensor. In a properly functioning satellite,
failures of this kind should be rare, but significant when they do
happen. Therefore we enabled the checking of this flag by
default."""),
    arcGISDisplayName=_(u'NAVWARN - Navigation quality is reduced'),
    arcGISCategory=_(u'l2_flags to check'))

AddArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkABSAER',
    typeMetadata=BooleanTypeMetadata(),
    description=_(
u"""If both this option and the option to use l2_flags to unmask
daytime images are enabled, then the ABSAER flag will be checked, and
if that flag is set, the pixel will remain masked.

This option is disabled by default.

`Here <http://oceancolor.gsfc.nasa.gov/VALIDATION/flags.html>`_, NASA
describes the ABSAER flag as "possible absorbing aerosol (disabled)".
Because NASA marked this flag as "disabled" in their documentation, we
have disabled checking of it by default."""),
    arcGISDisplayName=_(u'ABSAER - Possible absorbing aerosol'),
    arcGISCategory=_(u'l2_flags to check'))

AddArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkMAXAERITER',
    typeMetadata=BooleanTypeMetadata(),
    description=_(
u"""If both this option and the option to use l2_flags to unmask
daytime images are enabled, then the MAXAERITER flag will be checked, and
if that flag is set, the pixel will remain masked.

This option is enabled by default.

`Here <http://oceancolor.gsfc.nasa.gov/VALIDATION/flags.html>`_, NASA
describes the MAXAERITER flag as "aerosol iterations exceeded max". We
have not researched exactly what this problem means, but presume that
some anomalous or unexpected atmospheric conditions lead to a failure
in correcting for the possible presence of aerosols. To be cautious,
we assumed this could produce an error in the SST estimate, therefore
we enabled checking this flag by default."""),
    arcGISDisplayName=_(u'MAXAERITER - Aerosol iterations exceeded max'),
    arcGISCategory=_(u'l2_flags to check'))

AddArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkMODGLINT',
    typeMetadata=BooleanTypeMetadata(),
    description=_(
u"""If both this option and the option to use l2_flags to unmask
daytime images are enabled, then the MODGLINT flag will be checked, and
if that flag is set, the pixel will remain masked.

This option is disabled by default.

`Here <http://oceancolor.gsfc.nasa.gov/VALIDATION/flags.html>`_, NASA
describes the MODGLINT flag as "moderate sun glint contamination". As
we understand it, the longwave SST algorithm used for daytime SST
images is not affected by sun glint, so we disabled the checking of
this flag by default. For an image that illustrates this, see page 14
of `this presentation <http://oceancolor.gsfc.nasa.gov/DOCS/Presentations/seadas_umbc_franz.pdf>`_."""),
    arcGISDisplayName=_(u'MODGLINT - Moderate sun glint contamination'),
    arcGISCategory=_(u'l2_flags to check'))

AddArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkCHLWARN',
    typeMetadata=BooleanTypeMetadata(),
    description=_(
u"""If both this option and the option to use l2_flags to unmask
daytime images are enabled, then the CHLWARN flag will be checked, and
if that flag is set, the pixel will remain masked.

This option is disabled by default.

`Here <http://oceancolor.gsfc.nasa.gov/VALIDATION/flags.html>`_, NASA
describes the CHLWARN flag as "derived product quality is reduced".
Presumably the derived product is chlorophyll. As with the CHLFAIL
flag, we see no reason why SST estimation would be affected by
problems that are very specific to cholorophyll estimation, therefore
we disabled the checking of this flag by default."""),
    arcGISDisplayName=_(u'CHLWARN - Derived product quality is reduced'),
    arcGISCategory=_(u'l2_flags to check'))

AddArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkATMWARN',
    typeMetadata=BooleanTypeMetadata(),
    description=_(
u"""If both this option and the option to use l2_flags to unmask
daytime images are enabled, then the ATMWARN flag will be checked, and
if that flag is set, the pixel will remain masked.

This option is enabled by default.

`Here <http://oceancolor.gsfc.nasa.gov/VALIDATION/flags.html>`_, NASA
describes the ATMWARN flag as "atmospheric correction is suspect". We
did not research the circumstances under which this failure might
occur; the details might appear in
'this document <http://oceancolor.gsfc.nasa.gov/SeaWiFS/TECH_REPORTS/PLVol22.pdf>_.
As with the ATMFAIL flag, we assumed that if the atmospheric
correction correction was suspect, it was likely that the atmosphere
was in an anomalous or unexpected state, or that sufficient data were
not available to correct for the state of atmosphere, and that the SST
estimates would also have significant error. This may be a bad
assumption if the atmospheric correction algorithm is designed
explicitly for chlorophyll but we elected to be conservative and leave
the pixels masked."""),
    arcGISDisplayName=_(u'ATMWARN - Atmospheric correction is suspect'),
    arcGISCategory=_(u'l2_flags to check'))

AddArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkSEAICE',
    typeMetadata=BooleanTypeMetadata(),
    description=_(
u"""If both this option and the option to use l2_flags to unmask
daytime images are enabled, then the SEAICE flag will be checked, and
if that flag is set, the pixel will remain masked.

This option is enabled by default.

`Here <http://oceancolor.gsfc.nasa.gov/VALIDATION/flags.html>`_, NASA
describes the SEAICE flag as "possible sea ice contamination".
"""),
    arcGISDisplayName=_(u'SEAICE - Possible sea ice contamination'),
    arcGISCategory=_(u'l2_flags to check'))

AddArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkNAVFAIL',
    typeMetadata=BooleanTypeMetadata(),
    description=_(
u"""If both this option and the option to use l2_flags to unmask
daytime images are enabled, then the NAVFAIL flag will be checked, and
if that flag is set, the pixel will remain masked.

This option is enabled by default.

`Here <http://oceancolor.gsfc.nasa.gov/VALIDATION/flags.html>`_, NASA
describes the NAVFAIL flag as "bad navigation". In satellite remote
sensing, navigation refers to the problem of determining accurate and
precise geographic coordinates for each measurement taken by the
sensor. In a properly functioning satellite, failures of this kind
should be rare, but significant when they do happen. Therefore we
enabled the checking of this flag by default."""),
    arcGISDisplayName=_(u'NAVFAIL - Bad navigation'),
    arcGISCategory=_(u'l2_flags to check'))

AddArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkFILTER',
    typeMetadata=BooleanTypeMetadata(),
    description=_(
u"""If both this option and the option to use l2_flags to unmask
daytime images are enabled, then the FILTER flag will be checked, and
if that flag is set, the pixel will remain masked.

This option is disabled by default.

`Here <http://oceancolor.gsfc.nasa.gov/VALIDATION/flags.html>`_, NASA
describes the FILTER flag as "pixel rejected by user-defined filter".
We found no documentation for this flag. We assumed it was very
specific to chlorophyll processing and decided to ignore it for SST
processing, therefore we disabled the checking of this flag by
default."""),
    arcGISDisplayName=_(u'FILTER - Pixel rejected by user-defined filter'),
    arcGISCategory=_(u'l2_flags to check'))

AddArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkHIPOL',
    typeMetadata=BooleanTypeMetadata(),
    description=_(
u"""If both this option and the option to use l2_flags to unmask
daytime images are enabled, then the HIPOL flag will be checked, and
if that flag is set, the pixel will remain masked.

This option is disabled by default.

`Here <http://oceancolor.gsfc.nasa.gov/VALIDATION/flags.html>`_, NASA
describes the HIPOL flag as "high degree of polarization". We do not
know of any reason why a high degree of polarization would affect SST
the SST estimation, so we disabled this flag by default."""),
    arcGISDisplayName=_(u'HIPOL - High degree of polarization'),
    arcGISCategory=_(u'l2_flags to check'))

AddArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkPRODFAIL',
    typeMetadata=BooleanTypeMetadata(),
    description=_(
u"""If both this option and the option to use l2_flags to unmask
daytime images are enabled, then the PRODFAIL flag will be checked,
and if that flag is set, the pixel will remain masked.

This option is disabled by default.

`Here <http://oceancolor.gsfc.nasa.gov/VALIDATION/flags.html>`_, NASA
describes the PRODFAIL flag as "derived product failure". Presumably
the derived product is chlorophyll. As with the CHLFAIL flag, we see
no reason why SST estimation would be affected by problems that are
very specific to cholorophyll estimation, therefore we disabled the
checking of this flag by default."""),
    arcGISDisplayName=_(u'PRODFAIL - Derived product failure'),
    arcGISCategory=_(u'l2_flags to check'))

AddArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'overwriteExisting',
    typeMetadata=BooleanTypeMetadata(),
    description=_(
u"""If True, the output feature class will be overwritten, if it
exists. If False, a ValueError will be raised if the output feature
class exists."""),
    initializeToArcGISGeoprocessorVariable=u'OverwriteOutput')

# Public method: OceanColorLevel2Converter.L2LacSstHdfToArcGISRaster

AddMethodMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISRaster,
    shortDescription=_(u'Converts a NASA OceanColor MODIS L2 LAC SST file to an ArcGIS raster.'),
    longDescription=_(u'TODO: Write long description.'),
    isExposedToPythonCallers=True,
    isExposedByCOM=True,
    isExposedAsArcGISTool=True,
    arcGISDisplayName=_(u'Convert MODIS L2 LAC SST to ArcGIS Raster'),
    arcGISToolCategory=_(u'Data Products\\NASA GSFC OceanColor Group\\L2 Products'),
    dependencies=[ArcGISDependency(9, 1), ArcGISExtensionDependency(u'spatial'), PythonAggregatedModuleDependency('numpy')])

CopyArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'cls', OceanColorLevel2Converter.L2LacSstHdfToArcGISRaster, u'cls')
CopyArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'hdf', OceanColorLevel2Converter.L2LacSstHdfToArcGISRaster, u'hdf')

AddArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISRaster, u'raster',
    typeMetadata=ArcGISRasterTypeMetadata(mustBeDifferentThanArguments=[u'hdf'], deleteIfParameterIsTrue=u'overwriteExisting', createParentDirectories=True),
    description=_(
u"""TODO: Write description"""),
    direction=u'Output',
    arcGISDisplayName=_(u'Output raster'))

CopyArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'variable', OceanColorLevel2Converter.L2LacSstHdfToArcGISRaster, u'variable')
CopyArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'spatialExtent', OceanColorLevel2Converter.L2LacSstHdfToArcGISRaster, u'spatialExtent')
CopyArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'qualityLevel', OceanColorLevel2Converter.L2LacSstHdfToArcGISRaster, u'qualityLevel')
CopyArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'recoverBadPixels', OceanColorLevel2Converter.L2LacSstHdfToArcGISRaster, u'recoverBadPixels')
CopyArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkATMFAIL', OceanColorLevel2Converter.L2LacSstHdfToArcGISRaster, u'checkATMFAIL')
CopyArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkLAND', OceanColorLevel2Converter.L2LacSstHdfToArcGISRaster, u'checkLAND')
CopyArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkPRODWARN', OceanColorLevel2Converter.L2LacSstHdfToArcGISRaster, u'checkPRODWARN')
CopyArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkHIGLINT', OceanColorLevel2Converter.L2LacSstHdfToArcGISRaster, u'checkHIGLINT')
CopyArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkHILT', OceanColorLevel2Converter.L2LacSstHdfToArcGISRaster, u'checkHILT')
CopyArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkHISATZEN', OceanColorLevel2Converter.L2LacSstHdfToArcGISRaster, u'checkHISATZEN')
CopyArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkCOASTZ', OceanColorLevel2Converter.L2LacSstHdfToArcGISRaster, u'checkCOASTZ')
CopyArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkSTRAYLIGHT', OceanColorLevel2Converter.L2LacSstHdfToArcGISRaster, u'checkSTRAYLIGHT')
CopyArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkCLDICE', OceanColorLevel2Converter.L2LacSstHdfToArcGISRaster, u'checkCLDICE')
CopyArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkCOCCOLITH', OceanColorLevel2Converter.L2LacSstHdfToArcGISRaster, u'checkCOCCOLITH')
CopyArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkTURBIDW', OceanColorLevel2Converter.L2LacSstHdfToArcGISRaster, u'checkTURBIDW')
CopyArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkHISOLZEN', OceanColorLevel2Converter.L2LacSstHdfToArcGISRaster, u'checkHISOLZEN')
CopyArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkLOWLW', OceanColorLevel2Converter.L2LacSstHdfToArcGISRaster, u'checkLOWLW')
CopyArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkCHLFAIL', OceanColorLevel2Converter.L2LacSstHdfToArcGISRaster, u'checkCHLFAIL')
CopyArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkNAVWARN', OceanColorLevel2Converter.L2LacSstHdfToArcGISRaster, u'checkNAVWARN')
CopyArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkABSAER', OceanColorLevel2Converter.L2LacSstHdfToArcGISRaster, u'checkABSAER')
CopyArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkMAXAERITER', OceanColorLevel2Converter.L2LacSstHdfToArcGISRaster, u'checkMAXAERITER')
CopyArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkMODGLINT', OceanColorLevel2Converter.L2LacSstHdfToArcGISRaster, u'checkMODGLINT')
CopyArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkCHLWARN', OceanColorLevel2Converter.L2LacSstHdfToArcGISRaster, u'checkCHLWARN')
CopyArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkATMWARN', OceanColorLevel2Converter.L2LacSstHdfToArcGISRaster, u'checkATMWARN')
CopyArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkSEAICE', OceanColorLevel2Converter.L2LacSstHdfToArcGISRaster, u'checkSEAICE')
CopyArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkNAVFAIL', OceanColorLevel2Converter.L2LacSstHdfToArcGISRaster, u'checkNAVFAIL')
CopyArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkFILTER', OceanColorLevel2Converter.L2LacSstHdfToArcGISRaster, u'checkFILTER')
CopyArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkHIPOL', OceanColorLevel2Converter.L2LacSstHdfToArcGISRaster, u'checkHIPOL')
CopyArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISPoints, u'checkPRODFAIL', OceanColorLevel2Converter.L2LacSstHdfToArcGISRaster, u'checkPRODFAIL')

AddArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISRaster, u'templateRaster',
    typeMetadata=ArcGISRasterTypeMetadata(mustBeDifferentThanArguments=[u'hdf', u'raster'], canBeNone=True),
    description=_(
u"""TODO: Write description"""),
    arcGISDisplayName=_(u'Template raster'),
    arcGISCategory=_(u'Interpolation options'))

AddArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISRaster, u'power',
    typeMetadata=FloatTypeMetadata(mustBeGreaterThan=0.),
    description=_(
u"""TODO: Write description"""),
    arcGISDisplayName=_(u'Power'),
    arcGISCategory=_(u'Interpolation options'))

AddArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISRaster, u'searchRadius',
    typeMetadata=UnicodeStringTypeMetadata(allowedValues=[u'Fixed', u'Variable'], makeLowercase=True),
    description=_(
u"""TODO: Write description"""),
    arcGISDisplayName=_(u'Search radius'),
    arcGISCategory=_(u'Interpolation options'))

AddArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISRaster, u'distance',
    typeMetadata=FloatTypeMetadata(mustBeGreaterThan=0.),
    description=_(
u"""TODO: Write description"""),
    arcGISDisplayName=_(u'Search distance'),
    arcGISCategory=_(u'Interpolation options'))

AddArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISRaster, u'numPoints',
    typeMetadata=IntegerTypeMetadata(minValue=0),
    description=_(
u"""TODO: Write description"""),
    arcGISDisplayName=_(u'Number of points'),
    arcGISCategory=_(u'Interpolation options'))

AddArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISRaster, u'barrier',
    typeMetadata=ArcGISFeatureLayerTypeMetadata(canBeNone=True, mustExist=True, allowedShapeTypes=[u'Polyline']),
    description=_(
u"""TODO: Write description"""),
    arcGISDisplayName=_(u'Input barrier polyline features'),
    arcGISCategory=_(u'Interpolation options'))

AddArgumentMetadata(OceanColorLevel2Converter.L2LacSstHdfToArcGISRaster, u'overwriteExisting',
    typeMetadata=BooleanTypeMetadata(),
    description=_(
u"""If True, the output raster will be overwritten, if it exists. If
False, a ValueError will be raised if the output raster exists."""),
    initializeToArcGISGeoprocessorVariable=u'OverwriteOutput')

###############################################################################
# Names exported by this module
###############################################################################

__all__ = ['OceanColorLevel3SMIFileSearcher',
           'OceanColorLevel3SMITimeSeries',
           'OceanColorLevel2Converter']
