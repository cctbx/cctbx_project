from __future__ import division
from __future__ import print_function

class phil_validation:
  def __init__(self,param):

    self.param = param
    self.application_level_validation()

  def application_level_validation(self):
    pass

from xfel.merging.database.merging_database import manager
class application(manager):
  def __init__(self,params):
    self.params = params
    self.db = self.connection()
    self.cursor = self.db.cursor()
    self.insert_isoform()
    self.insert_hkl()
    self.db.commit()

  def connection(self):
    print()
    if self.params.db.password is None:
      import getpass
      password = getpass.getpass()
    else:
      password = self.params.db.password

    try:
      import MySQLdb
      db = MySQLdb.connect(passwd=password,
                           user = self.params.db.user,
                           host = self.params.db.host,
                           db = self.params.db.name,compress=False)
      cursor = db.cursor()
      cursor.execute("use %s;"%self.params.db.name)

      return db
    except Exception as e:
      print(e)
      raise RuntimeError("Couldn't connect to mysql database")

  def insert_isoform(self, **kwargs):

    kwargs["name"]=self.params.isoform.name
    cell = self.params.isoform.cell.parameters()
    kwargs["cell_a"]=cell[0]
    kwargs["cell_b"]=cell[1]
    kwargs["cell_c"]=cell[2]
    kwargs["cell_alpha"]=cell[3]
    kwargs["cell_beta"]=cell[4]
    kwargs["cell_gamma"]=cell[5]
    kwargs["lookup_symbol"]=self.params.isoform.lookup_symbol

    (sql, parameters) = self._insert(
      table='`%s_isoforms`' % self.params.experiment_tag,
      **kwargs)

    self.cursor.execute(sql, parameters[0])
    self._lastrowid = self.cursor.lastrowid
    print("isoformID=",self._lastrowid)

  def insert_hkl(self):
    kwargs = {}
    cell = self.params.isoform.cell
    from cctbx.crystal import symmetry

    cs = symmetry(unit_cell = cell,space_group_symbol=self.params.isoform.lookup_symbol)
    mset = cs.build_miller_set(anomalous_flag=False, d_min=self.params.isoform.resolution_limit)

    indices = mset.indices()

    from six.moves import cStringIO as StringIO
    query = StringIO()
    query.write("INSERT INTO `%s_hkls` (h,k,l,isoforms_isoform_id) VALUES "%self.params.experiment_tag)
    firstcomma = ""
    for item in indices:
      query.write(firstcomma); firstcomma=","
      query.write("('%d','%d','%d','%d')"%(item[0],item[1],item[2],self._lastrowid))
    self.cursor.execute( query.getvalue() )

    print("Inserted %d HKLs"%(len(indices)))
