from __future__ import division
from __future__ import print_function
from cctbx.array_family import flex

mysql_master_phil = """
backend = FS *MySQL SQLite Flex
  .type = choice
  .help = "Back end database; FS for flat-file ASCII data storage,
           MySQL and SQLite for the respective proper database
           backends. Flex gives in-memory flex arrays instead of disk-based storage,
           which are ultimately written as pickle files at the final join()"
mysql {
  # MySQL database v5.1 data store.
  # mysql -u root -p # Steps to be taken by the database administrator
  # CREATE DATABASE database;
  # GRANT ALL ON database.* to 'user' IDENTIFIED BY 'passwd';
  # SET GLOBAL max_allowed_packet=512*1024*1024;
  # installation of MySQLdb, download from http://sourceforge.net/projects/mysql-python
  # install into the cctbx python with libtbx.python setup.py install
  # Maintenance and cleanup by user with mysql -u user -p
  # SHOW TABLES FROM database;
  # DROP TABLE *;
  host = localhost
    .type = str
    .help = persistent data tables to MySQL database using mysql-server on this host
    .help = concurrent client connections OK, can use nproc > 1
  port = 3306
    .type = int
    .help = port number for connecting on this host
  runtag = None
    .type = str
    .help = "Identifier is for this run, signifies database to use.
             Dumps old data when applicable and writes new tables."
  user = None
    .type = str
    .help = mysql username provided by the database administrator
  passwd = None
    .type = str
    .help = mysql password provided by the database administrator
  database = None
    .type = str
    .help = mysql user's working database name, provided by the database administrator
}
"""
class manager_base (object):
  def insert_frame_legacy(self,result,wavelength,corr,slope,offset,data):
    from scitbx import matrix
    """Legacy compatibility with cxi.merge; insert frame-data to backend.
    XXX needs to be backported to SQLite backend (use this base class)
    result: an unpickled dictionary from an integration pickle
    wavelength, beam_x, beam_y, distance: parameters from the model
    data: an instance of a "frame_data" container class
    postx: an instance of the legacy_cxi_merge_postrefinement results, or None
    """
    have_sa_params = ( type(result.get("sa_parameters")[0]) == type(dict()) )

    kwargs = {'wavelength': wavelength,
              'beam_x': result['xbeam'],
              'beam_y': result['ybeam'],
              'distance': result['distance'],
              'c_c': corr,
              'slope': slope,
              'offset': offset,
              'unique_file_name': data.file_name}
    if have_sa_params:
      sa_parameters = result['sa_parameters'][0]
      res_ori_direct = sa_parameters['reserve_orientation'].direct_matrix().elems

      kwargs['res_ori_1'] = res_ori_direct[0]
      kwargs['res_ori_2'] = res_ori_direct[1]
      kwargs['res_ori_3'] = res_ori_direct[2]
      kwargs['res_ori_4'] = res_ori_direct[3]
      kwargs['res_ori_5'] = res_ori_direct[4]
      kwargs['res_ori_6'] = res_ori_direct[5]
      kwargs['res_ori_7'] = res_ori_direct[6]
      kwargs['res_ori_8'] = res_ori_direct[7]
      kwargs['res_ori_9'] = res_ori_direct[8]

      kwargs['rotation100_rad'] = sa_parameters.rotation100_rad
      kwargs['rotation010_rad'] = sa_parameters.rotation010_rad
      kwargs['rotation001_rad'] = sa_parameters.rotation001_rad

      kwargs['half_mosaicity_deg'] = sa_parameters.half_mosaicity_deg
      kwargs['wave_HE_ang'] = sa_parameters.wave_HE_ang
      kwargs['wave_LE_ang'] = sa_parameters.wave_LE_ang
      kwargs['domain_size_ang'] = sa_parameters.domain_size_ang

    else:
      res_ori_direct = matrix.sqr(
        data.indexed_cell.orthogonalization_matrix()).transpose().elems

      kwargs['res_ori_1'] = res_ori_direct[0]
      kwargs['res_ori_2'] = res_ori_direct[1]
      kwargs['res_ori_3'] = res_ori_direct[2]
      kwargs['res_ori_4'] = res_ori_direct[3]
      kwargs['res_ori_5'] = res_ori_direct[4]
      kwargs['res_ori_6'] = res_ori_direct[5]
      kwargs['res_ori_7'] = res_ori_direct[6]
      kwargs['res_ori_8'] = res_ori_direct[7]
      kwargs['res_ori_9'] = res_ori_direct[8]
      if self.params.scaling.report_ML:
        kwargs['half_mosaicity_deg'] = result["ML_half_mosaicity_deg"][0]
        kwargs['domain_size_ang'] = result["ML_domain_size_ang"][0]
      else:
        kwargs['half_mosaicity_deg'] =float("NaN")
        kwargs['domain_size_ang'] =float("NaN")
    return self.insert_frame(**kwargs)

class manager (manager_base):
  def __init__(self,params):
    self.params = params

  def connection(self):
    try:
      import MySQLdb
      db = MySQLdb.connect(passwd=self.params.mysql.passwd,
                           user = self.params.mysql.user,
                           host = self.params.mysql.host,
                           port = self.params.mysql.port,
                           db = self.params.mysql.database,compress=False)
      cursor = db.cursor()
      cursor.execute("use %s;"%self.params.mysql.database)
      db.commit()

      return db
    except Exception:
      raise RuntimeError("Couldn't connect to mysql database")


  def initialize_db(self, indices):
    db = self.connection()
    print("testing for tables")
    cursor = db.cursor()
    cursor.execute("SHOW TABLES from %s;"%self.params.mysql.database)
    all_tables = cursor.fetchall()

    # Beware of SQL injection vulnerability (here and elsewhere).
    new_tables = self.merging_schema_tables(self.params.mysql.runtag)
    for table in new_tables:
      cursor.execute("DROP TABLE IF EXISTS %s;"%table[0])
      cursor.execute("CREATE TABLE %s "%table[0]+table[1].replace("\n"," ")+" ;")
    from six.moves import cStringIO as StringIO
    query = StringIO()
    query.write("INSERT INTO `%s_miller` (h,k,l) VALUES "%self.params.mysql.runtag)
    firstcomma = ""
    for item in indices:
      query.write(firstcomma); firstcomma=","
      query.write("('%d','%d','%d')"%(item[0],item[1],item[2]))
    query.write(" ;")
    cursor.execute( query.getvalue() )
    db.commit()

  def _insert(self, table, **kwargs):
    """The _insert() function generates the SQL command and parameter
    argument for the _execute() function.
    """

    # Note MySQL uses "%s", whereas SQLite uses "?".
    sql = ("INSERT INTO %s (" % table) \
          + ", ".join(kwargs.keys()) + ") VALUES (" \
          + ", ".join(["%s"] * len(kwargs.keys())) + ")"

    # If there are more than one rows to insert, "unpack" the keyword
    # argument iterables and zip them up.  This effectively rearranges
    # a list of columns into a list of rows.
    try:
      parameters = zip(*kwargs.values())
    except TypeError:
      parameters = [kwargs.values()]

    return (sql, parameters)


  def insert_frame(self, **kwargs):
    db = self.connection()
    cursor = db.cursor()

    (sql, parameters) = self._insert(
      table='`%s_frame`' % self.params.mysql.runtag,
      **kwargs)

    cursor.execute(sql, parameters[0])
    db.commit()

    # Entry in the observation table is zero-based.
    return cursor.lastrowid - 1


  def insert_observation(self, **kwargs):
    return
    db = self.connection()
    cursor = db.cursor()

    # For MySQLdb executemany() is six times slower than a single big
    # execute() unless the "values" keyword is given in lowercase
    # (http://sourceforge.net/p/mysql-python/bugs/305).
    #
    # See also merging_database_sqlite3._insert()
    query = ("INSERT INTO `%s_observation` (" % self.params.mysql.runtag) \
            + ", ".join(kwargs.keys()) + ") values (" \
            + ", ".join(["%s"] * len(kwargs.keys())) + ")"
    try:
      parameters = zip(*kwargs.values())
    except TypeError:
      parameters = [kwargs.values()]
    cursor.executemany(query, parameters)
    db.commit()


  def join(self):
    pass


  def read_indices(self):
    db = self.connection()
    cursor = db.cursor()
    from cctbx.array_family import flex
    millers = dict(merged_asu_hkl=flex.miller_index())
    cursor.execute("SELECT h,k,l FROM `%s_miller` ORDER BY hkl_id_1_base"%self.params.mysql.runtag)
    for item in cursor.fetchall():
      millers["merged_asu_hkl"].append((item[0],item[1],item[2]))
    return millers

  def read_observations(self):
    db = self.connection()
    cursor = db.cursor()
    cursor.execute("SELECT hkl_id_0_base,i,sigi,frame_id_0_base,original_h,original_k,original_l FROM `%s_observation`"%self.params.mysql.runtag)
    ALL = cursor.fetchall()

    return dict(hkl_id = flex.int([a[0] for a in ALL]), #as MySQL indices are 1-based
               i = flex.double([a[1] for a in ALL]),
               sigi = flex.double([a[2] for a in ALL]),
               frame_id = flex.int([a[3] for a in ALL]),
               original_h = flex.int([a[4] for a in ALL]),
               original_k = flex.int([a[5] for a in ALL]),
               original_l = flex.int([a[6] for a in ALL]),
               )

  def read_frames(self):
    from xfel.cxi.util import is_odd_numbered
    db = self.connection()
    cursor = db.cursor()
    cursor.execute("""SELECT
    frame_id_1_base,wavelength,c_c,slope,offset,res_ori_1,res_ori_2,res_ori_3,
    res_ori_4,res_ori_5,res_ori_6,res_ori_7,res_ori_8,res_ori_9,
    unique_file_name
    FROM `%s_frame`"""%self.params.mysql.runtag)
    ALL = cursor.fetchall()
    from cctbx.crystal_orientation import crystal_orientation
    orientations = [crystal_orientation(
     (a[5],a[6],a[7],a[8],a[9],a[10],a[11],a[12],a[13]),False) for a in ALL]
    return dict( frame_id = flex.int( [a[0]-1 for a in ALL] ),
               wavelength = flex.double( [a[1] for a in ALL] ),
                       cc = flex.double( [a[2] for a in ALL] ),
                    slope = flex.double( [a[3] for a in ALL] ),
                   offset = flex.double( [a[4] for a in ALL] ),
             odd_numbered = flex.bool( [is_odd_numbered(a[14], use_hash = self.params.hash_filenames) for a in ALL] ),
              orientation = orientations,
                unit_cell = [CO.unit_cell() for CO in orientations],
         unique_file_name = [a[14] for a in ALL] )

  def merging_schema_tables(self,runtag):
    return [("`"+runtag+"_observation`","""
            (
              hkl_id_0_base INT,
              i DOUBLE NOT NULL,
              sigi DOUBLE NOT NULL,
              detector_x DOUBLE NOT NULL,
              detector_y DOUBLE NOT NULL,
              frame_id_0_base INT,
              overload_flag INTEGER,
              original_h INT NOT NULL,
              original_k INT NOT NULL,
              original_l INT NOT NULL
            )
            """),
            ("`"+runtag+"_frame`","""
            (
              frame_id_1_base INT UNSIGNED AUTO_INCREMENT NOT NULL PRIMARY KEY,
              wavelength DOUBLE NOT NULL,
              beam_x DOUBLE NOT NULL,
              beam_y DOUBLE NOT NULL,
              distance DOUBLE NOT NULL,
              c_c DOUBLE NOT NULL,
              slope DOUBLE NOT NULL,
              offset DOUBLE NOT NULL,
              res_ori_1 DOUBLE NOT NULL,
              res_ori_2 DOUBLE NOT NULL,
              res_ori_3 DOUBLE NOT NULL,
              res_ori_4 DOUBLE NOT NULL,
              res_ori_5 DOUBLE NOT NULL,
              res_ori_6 DOUBLE NOT NULL,
              res_ori_7 DOUBLE NOT NULL,
              res_ori_8 DOUBLE NOT NULL,
              res_ori_9 DOUBLE NOT NULL,
              rotation100_rad DOUBLE,
              rotation010_rad DOUBLE,
              rotation001_rad DOUBLE,
              half_mosaicity_deg DOUBLE,
              wave_HE_ang DOUBLE,
              wave_LE_ang DOUBLE,
              domain_size_ang DOUBLE,
              unique_file_name MEDIUMTEXT
              ) AUTO_INCREMENT = 1
            """
            ),
            ("`"+runtag+"_miller`","""(
              hkl_id_1_base INT AUTO_INCREMENT PRIMARY KEY,
              h INT NOT NULL,
              k INT NOT NULL,
              l INT NOT NULL
              ) AUTO_INCREMENT = 1
            """
            ),
              ]
  def positional_refinement_schema_tables(self,runtag):
    return [("`"+runtag+"_spotfinder`","""
            (
              frame_id INT, itile INT,
              beam1x DOUBLE NOT NULL,
              beam1y DOUBLE NOT NULL,
              beamrx DOUBLE NOT NULL,
              beamry DOUBLE NOT NULL,
              spotfx DOUBLE NOT NULL,
              spotfy DOUBLE NOT NULL,
              spotcx DOUBLE NOT NULL,
              spotcy DOUBLE NOT NULL,
              h INT NOT NULL,
              k INT NOT NULL,
              l INT NOT NULL,
              radialpx DOUBLE NOT NULL DEFAULT 0.0,
              azimutpx DOUBLE NOT NULL DEFAULT 0.0
            )
            """),
              ]
