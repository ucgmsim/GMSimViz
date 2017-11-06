import datetime
class shakemap_grid:
    def __init__(self, filename):
        self.fp = open(filename, 'w')

    def write_shakemap_grid_header(self, event_id, event_type, mag, dep, hlat, hlon, origin_time, run_name, x_min, x_max, y_min, y_max, grd_nx, grd_ny):
        self.fp.write('<?xml version="1.0" \
    encoding="UTF-8" standalone="yes"?>\n')
        self.fp.write('<shakemap_grid \
    xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" \
    xmlns="http://earthquake.usgs.gov/eqcenter/shakemap" \
    xsi:schemaLocation="http://earthquake.usgs.gov \
    http://earthquake.usgs.gov/eqcenter/shakemap/xml/schemas/shakemap.xsd" \
    event_id="%s" shakemap_id="%s" \
    shakemap_version="1" code_version="1" \
    process_timestamp="%s" \
    shakemap_originator="nz" \
    map_status="RELEASED" \
    shakemap_event_type="%s">\n' % (event_id, event_id, \
                                    datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ"), event_type))
        self.fp.write('<event event_id="%s" magnitude="%s" \
    depth="%s" lat="%s" lon="%s" event_timestamp="%s" \
    event_network="nz" event_description="%s" />\n' % (event_id, mag, \
                                                       dep, hlat, hlon, origin_time, run_name.split('_')[0]))
        self.fp.write('<grid_specification lon_min="%s" lat_min="%s" \
    lon_max="%s" lat_max="%s" nominal_lon_spacing="%s" \
    nominal_lat_spacing="%s" nlon="%d" nlat="%d" />\n' % (x_min, y_min, \
                                                          x_max, y_max, abs(x_max - x_min) / grd_nx, \
                                                          abs(y_max - y_min) / grd_ny, grd_nx, grd_ny))
        self.fp.write('<event_specific_uncertainty name="pgv" \
    value="-1" numsta="0" />\n')
        self.fp.write('<event_specific_uncertainty name="mi" \
    value="-1" numsta="0" />\n')
        self.fp.write('<grid_field index="1" name="LON" units="dd" />\n')
        self.fp.write('<grid_field index="2" name="LAT" units="dd" />\n')
        self.fp.write('<grid_field index="3" name="PGV" units="cms" />\n')
        self.fp.write('<grid_field index="4" name="MMI" units="intensity" />\n')
        self.fp.write('<grid_data>\n')

    def write(self, string):
        self.fp.write(string)

    def write_shakemap_grid_footer():
        self.fp.write('</grid_data>\n</shakemap_grid>\n')
        self.fp.close()