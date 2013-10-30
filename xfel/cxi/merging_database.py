from __future__ import division
from cctbx.array_family import flex

mysql_master_phil = """
backend = FS *MySQL SQLite
  .type = choice
  .help = "Back end database; FS for flat-file ASCII data storage,
           MySQL and SQLite for the respective proper database
           backends."
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

class manager:
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

      return db
    except Exception:
      raise RuntimeError("Couldn't connect to mysql database")


  def initialize_db(self, indices):
    db = self.connection()
    print "testing for tables"
    cursor = db.cursor()
    cursor.execute("SHOW TABLES from %s;"%self.params.mysql.database)
    all_tables = cursor.fetchall()

    # Beware of SQL injection vulnerability (here and elsewhere).
    new_tables = self.merging_schema_tables(self.params.mysql.runtag)
    for table in new_tables:
      cursor.execute("DROP TABLE IF EXISTS %s;"%table[0])
      cursor.execute("CREATE TABLE %s "%table[0]+table[1].replace("\n"," ")+" ;")
    import cStringIO
    query = cStringIO.StringIO()
    query.write("INSERT INTO %s_miller (h,k,l) VALUES "%self.params.mysql.runtag)
    firstcomma = ""
    for item in indices:
      query.write(firstcomma); firstcomma=","
      query.write("('%d','%d','%d')"%(item[0],item[1],item[2]))
    cursor.execute( query.getvalue() )


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
      table='%s_frame' % self.params.mysql.runtag,
      **kwargs)

    cursor.execute(sql, parameters[0])

    # Entry in the observation table is zero-based.
    return cursor.lastrowid - 1


  def insert_observation(self, **kwargs):
    db = self.connection()
    cursor = db.cursor()

    # For MySQLdb, executemany() is six times slower than the single
    # big command below.
    import cStringIO
    query = cStringIO.StringIO()
    query.write("""INSERT INTO %s_observation
    (hkl_id_0_base,i,sigi,detector_x,detector_y,frame_id_0_base,overload_flag,original_h,original_k,original_l)
    VALUES """%self.params.mysql.runtag)
    firstcomma = ""
    for i in range(len(kwargs['hkl_id_0_base'])):
      query.write(firstcomma); firstcomma=","
      query.write("('%7d','%18.8f','%18.8f','%8.2f','%8.2f','%7d','%d','%d','%d','%d')"%(
        kwargs['hkl_id_0_base'][i],
        kwargs['i'][i],
        kwargs['sigi'][i],
        kwargs['detector_x'][i],
        kwargs['detector_y'][i],
        kwargs['frame_id_0_base'][i],
        kwargs['overload_flag'][i],
        kwargs['original_h'][i],
        kwargs['original_k'][i],
        kwargs['original_l'][i]))
    cursor.execute(query.getvalue())


  def join(self):
    pass


  def read_indices(self):
    db = self.connection()
    cursor = db.cursor()
    from cctbx.array_family import flex
    millers = dict(merged_asu_hkl=flex.miller_index())
    cursor.execute("SELECT h,k,l FROM %s_miller ORDER BY hkl_id_1_base"%self.params.mysql.runtag)
    for item in cursor.fetchall():
      millers["merged_asu_hkl"].append((item[0],item[1],item[2]))
    return millers

  def read_observations(self):
    db = self.connection()
    cursor = db.cursor()
    cursor.execute("SELECT hkl_id_0_base,i,sigi,frame_id_0_base,original_h,original_k,original_l FROM %s_observation"%self.params.mysql.runtag)
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
    FROM %s_frame"""%self.params.mysql.runtag)
    ALL = cursor.fetchall()
    from cctbx.crystal_orientation import crystal_orientation
    orientations = [crystal_orientation(
     (a[5],a[6],a[7],a[8],a[9],a[10],a[11],a[12],a[13]),False) for a in ALL]
    return dict( frame_id = flex.int( [a[0]-1 for a in ALL] ),
               wavelength = flex.double( [a[1] for a in ALL] ),
                       cc = flex.double( [a[2] for a in ALL] ),
                    slope = flex.double( [a[3] for a in ALL] ),
                   offset = flex.double( [a[4] for a in ALL] ),
             odd_numbered = flex.bool( [is_odd_numbered(a[14]) for a in ALL] ),
              orientation = orientations,
                unit_cell = [CO.unit_cell() for CO in orientations] )

  def merging_schema_tables(self,runtag):
    return [(runtag+"_observation","""
            (
              hkl_id_0_base INT,
              i DOUBLE(18,8) NOT NULL,
              sigi DOUBLE(18,8) NOT NULL,
              detector_x DOUBLE(8,2) NOT NULL,
              detector_y DOUBLE(8,2) NOT NULL,
              frame_id_0_base INT,
              overload_flag INTEGER,
              original_h INT NOT NULL,
              original_k INT NOT NULL,
              original_l INT NOT NULL
            )
            """),
            (runtag+"_frame","""
            (
              frame_id_1_base INT UNSIGNED AUTO_INCREMENT NOT NULL PRIMARY KEY,
              wavelength DOUBLE(14,8) NOT NULL,
              beam_x DOUBLE(14,8) NOT NULL,
              beam_y DOUBLE(14,8) NOT NULL,
              distance DOUBLE(14,8) NOT NULL,
              c_c DOUBLE(10,7) NOT NULL,
              slope DOUBLE(11,8) NOT NULL,
              offset DOUBLE(10,2) NOT NULL,
              res_ori_1 DOUBLE(14,8) NOT NULL,
              res_ori_2 DOUBLE(14,8) NOT NULL,
              res_ori_3 DOUBLE(14,8) NOT NULL,
              res_ori_4 DOUBLE(14,8) NOT NULL,
              res_ori_5 DOUBLE(14,8) NOT NULL,
              res_ori_6 DOUBLE(14,8) NOT NULL,
              res_ori_7 DOUBLE(14,8) NOT NULL,
              res_ori_8 DOUBLE(14,8) NOT NULL,
              res_ori_9 DOUBLE(14,8) NOT NULL,
              rotation100_rad DOUBLE(10,7),
              rotation010_rad DOUBLE(10,7),
              rotation001_rad DOUBLE(10,7),
              half_mosaicity_deg DOUBLE(10,7),
              wave_HE_ang DOUBLE(14,8),
              wave_LE_ang DOUBLE(14,8),
              domain_size_ang DOUBLE(10,2),
              unique_file_name MEDIUMTEXT
              ) AUTO_INCREMENT = 1
            """
            ),
            (runtag+"_miller","""(
              hkl_id_1_base INT AUTO_INCREMENT PRIMARY KEY,
              h INT NOT NULL,
              k INT NOT NULL,
              l INT NOT NULL
              ) AUTO_INCREMENT = 1
            """
            ),
              ]
  def positional_refinement_schema_tables(self,runtag):
    return [(runtag+"_spotfinder","""
            (
              frame_id INT, itile INT,
              beam1x DOUBLE(10,2) NOT NULL,
              beam1y DOUBLE(10,2) NOT NULL,
              beamrx DOUBLE(10,2) NOT NULL,
              beamry DOUBLE(10,2) NOT NULL,
              spotfx DOUBLE(10,2) NOT NULL,
              spotfy DOUBLE(10,2) NOT NULL,
              spotcx DOUBLE(10,2) NOT NULL,
              spotcy DOUBLE(10,2) NOT NULL,
              h INT NOT NULL,
              k INT NOT NULL,
              l INT NOT NULL,
              radialpx DOUBLE(6,3) NOT NULL DEFAULT 0.0,
              azimutpx DOUBLE(6,3) NOT NULL DEFAULT 0.0
            )
            """),
              ]
