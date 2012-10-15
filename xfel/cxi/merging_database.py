from __future__ import division
from cctbx.array_family import flex
class manager:
  def __init__(self,params):
    self.params = params
  def use_mysql(self):
    return self.params.mysql.runtag is not None

  def connection(self):
    try:
      assert self.use_mysql()
      import MySQLdb
      db = MySQLdb.connect(passwd=self.params.mysql.passwd,
                           user = self.params.mysql.user,
                           host = self.params.mysql.host,
                           db = self.params.mysql.database,compress=False)
      cursor = db.cursor()
      cursor.execute("use %s;"%self.params.mysql.database)

      return db
    except Exception:
      print "Couldn't connect to mysql database"
      self.params.mysql.runtag = None

  def initialize_tag(self):
    db = self.connection()
    assert self.use_mysql()
    print "testing for tables"
    cursor = db.cursor()
    cursor.execute("SHOW TABLES from %s;"%self.params.mysql.database)
    all_tables = cursor.fetchall()

    new_tables = self.merging_schema_tables(self.params.mysql.runtag)
    for table in new_tables:
      cursor.execute("DROP TABLE IF EXISTS %s;"%table[0])
      cursor.execute("CREATE TABLE %s "%table[0]+table[1].replace("\n"," ")+" ;")

  def fill_indices(self,idxs):
    db = self.connection()
    assert self.use_mysql()
    cursor = db.cursor()
    import cStringIO
    query = cStringIO.StringIO()
    query.write("INSERT INTO %s_miller (h,k,l) VALUES "%self.params.mysql.runtag)
    firstcomma = ""
    for item in idxs:
      query.write(firstcomma); firstcomma=","
      query.write("('%d','%d','%d')"%(item[0],item[1],item[2]))
    cursor.execute( query.getvalue() )

  def read_indices(self):
    db = self.connection()
    assert self.use_mysql()
    cursor = db.cursor()
    from cctbx.array_family import flex
    millers = dict(merged_asu_hkl=flex.miller_index())
    cursor.execute("SELECT h,k,l FROM %s_miller ORDER BY hkl_id_1_base"%self.params.mysql.runtag)
    for item in cursor.fetchall():
      millers["merged_asu_hkl"].append((item[0],item[1],item[2]))
    return millers

  def read_observations(self):
    db = self.connection()
    assert self.use_mysql()
    cursor = db.cursor()
    cursor.execute("SELECT hkl_id_0_base,i,sigi,frame_id_0_base FROM %s_observation"%self.params.mysql.runtag)
    ALL = cursor.fetchall()

    return dict(hkl_id = flex.int([a[0] for a in ALL]), #as MySQL indices are 1-based
               i = flex.double([a[1] for a in ALL]),
               sigi = flex.double([a[2] for a in ALL]),
               frame_id = flex.int([a[3] for a in ALL]))

  def read_frames(self):
    from xfel.cxi.util import is_odd_numbered
    db = self.connection()
    assert self.use_mysql()
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
              i DOUBLE(14,8) NOT NULL,
              sigi DOUBLE(14,8) NOT NULL,
              detector_x DOUBLE(8,2) NOT NULL,
              detector_y DOUBLE(8,2) NOT NULL,
              frame_id_0_base INT,
              overload_flag ENUM('T','F')
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
