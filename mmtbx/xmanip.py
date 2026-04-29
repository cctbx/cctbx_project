"""Manipulation of xray data objects (coordinates and reflection files
"""

from __future__ import absolute_import, division, print_function
from mmtbx import xmanip_tasks
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
from iotbx import crystal_symmetry_from_any
import iotbx.phil
import iotbx.pdb
from cctbx import crystal
from cctbx import sgtbx
from cctbx import adptbx
from cctbx.array_family import flex
from scitbx.math import matrix
from libtbx.utils import Sorry, multi_out
from six.moves import cStringIO as StringIO
import sys, os
from six.moves import zip
from six.moves import range


class quick_rt_mx(object):
  from scitbx import matrix
  def __init__(self,
               r = [1,0,0,0,1,0,0,0,1] ,
               t = (0,0,0) ):
    self.rm = matrix.sqr(r)
    self.tv = flex.double(t)
  def r(self):
    return self.rm
  def t(self):
    return self.tv
  def inverse(self):
    tmp = quick_rt_mx(r=self.rm.inverse(), t=(-self.tv[0], -self.tv[1], -self.tv[2]) )
    return tmp
  def show(self,out=None):
    if out is None:
      out=sys.stdout
    print(file=out)
    print("R :", file=out)
    print(self.r().mathematica_form( one_row_per_line=True, format="%6.3f"), file=out)
    print("T :", file=out)
    for item in self.t():
      print("%6.3f"%(item), end=' ', file=out)
    print(file=out)
    print(file=out)

def construct_output_labels(labels, label_appendix, out=None ):
  if out is None:
    out = sys.stdout

  new_label_root = []
  standard_prefix = ["F","I","SIG", "HLA", "HL"]
  standard_postfix = ["(+)","(-)","PLUS", "MINUS","+","-","MINU" ]

  for label, app in zip(labels,label_appendix):
    if app is None:
      app=""
    #for each label, check if the prefix is present
    tmp_root = str(label[0])
    for post in standard_postfix:
      if tmp_root.endswith( post ):
        tmp_root = tmp_root[:len(post)-1]
        break

    for pre in standard_prefix:
      if tmp_root.startswith( pre ):
        tmp_root = tmp_root[len(pre)-1:]
        break

    if len(tmp_root)==0:
      tmp_root = label[0]

    if tmp_root[ len(tmp_root)-1 ]=="_":
      tmp_root = tmp_root[: len(tmp_root)-1 ]

    new_label_root.append( tmp_root+app )

  return new_label_root


def read_data(file_name, labels, xs, log=None):
  if log is None:
    log=sys.stdout
  if not os.path.isfile(file_name):
    raise Sorry("No such file: >%s<"%(file_name) )
  reflection_file = reflection_file_reader.any_reflection_file(
    file_name= file_name)
  miller_arrays = reflection_file.as_miller_arrays(crystal_symmetry=xs)
  label_table = reflection_file_utils.label_table(miller_arrays)
  return label_table.select_array(
    label=labels,
    command_line_switch="xray_data.labels",
    f=log)

def chain_name_modifier(chain_id, increment):
  ids=["A","B","C","D","E","F","G","H","I","J","K",
       "L","M","N","O","P","Q","R","S","T","U","V",
       "W","X","Y","Z","0","1","2","3","4","5","6",
       "7","8","9"]
  result=chain_id
  if chain_id in ids:
    new_index = (ids.index( chain_id ) + increment)%len(ids)
    result = ids[new_index]
  return result

def write_as_pdb_file(input_xray_structure = None,
                      input_crystal_symmetry = None,
                      input_pdb = None,
                      out = None,
                      chain_id_increment=5,
                      additional_remark=None,
                      print_cryst_and_scale=True):
  assert chain_id_increment is not None
  if out is None: out = sys.stdout

  xs = input_crystal_symmetry
  if xs is None: xs = input_xray_structure

  if additional_remark is not None:
    print("REMARK %s" % additional_remark, file=out)
  if print_cryst_and_scale:
    print("REMARK Space group: %s" % str(xs.space_group_info()), file=out)
  else:
    xs = None

  pdb_atoms = input_pdb.atoms()
  pdb_atoms.set_xyz(
    new_xyz=input_xray_structure.sites_cart())
  pdb_atoms.set_b(
    new_b=input_xray_structure.scatterers().extract_u_iso()
         *adptbx.u_as_b_factor)

  hierarchy = input_pdb.construct_hierarchy()
  for chain in hierarchy.chains():
    chain.id = chain_name_modifier(chain.id, chain_id_increment)

  print(hierarchy.as_pdb_string(crystal_symmetry=xs), file=out)

master_params = iotbx.phil.parse("""\
xmanip{
  input{
    unit_cell=None
    .type=unit_cell
    .help="Unit cell parameters"
    space_group=None
    .type=space_group
    .help="space group"
    xray_data
    .multiple=True
    .help="Scope defining xray data. Multiple scopes are allowed"
    {
      file_name=None
      .type=path
      .help="file name"
      labels=None
      .type=str
      .help="A unique label or unique substring of a label"
      label_appendix=None
      .type=str
      .help="Label appendix for output mtz file"
      name = None
      .type = str
      .help="An identifier of this particular miller array"
      write_out=None
      .type=bool
      .help="Determines if this data is written to the output file"

    }
    model
    .help="A model associated with the miller arrays. Only one model can be defined."
    {
      file_name=None
      .type=path
      .help="A model file"
    }
  }
  parameters{
    action = *reindex manipulate_pdb manipulate_miller
    .type=choice
    .help="Defines which action will be carried out."
    reindex
    .help="Reindexing parameters. Acts on coordinates and miller arrays."
    {
      standard_laws = niggli *reference_setting primitive_setting invert user_supplied
      .type=choice
      .help="Choices of reindexing operators. Will be applied on structure and miller arrays."
      user_supplied_law='h,k,l'
      .type=str
      .help="User supplied operator."
    }
    manipulate_miller
    .help="Acts on a single miller array or a set of miller arrays."
    {
      include scope mmtbx.xmanip_tasks.master_params
    }
    manipulate_pdb
    .help="Manipulate elements of a pdb file"
    {
      task = set_b apply_operator *None
      .type=choice
      .help="How to manipulate a pdb file"
      set_b{
        b_iso = 30
        .type=float
        .help="new B value for all atoms"
      }
      apply_operator{
        standard_operators = *user_supplied_operator user_supplied_cartesian_rotation_matrix
        .type=choice
        .help="Possible operators"
        user_supplied_operator = "x,y,z"
        .type=str
        .help="Actualy operator in x,y,z notation"
        user_supplied_cartesian_rotation_matrix
        .help="Rotation,translation matrix in cartesian frame"
        {
          r = None
          .type=str
          .help="Rotational part of operator"
          t = None
          .type = str
          .help="Translational part of operator"
        }
        invert = False
        .type = bool
        .help = "Invert operator given above before applying on coordinates"
        concatenate_model=False
        .type=bool
        .help="Determines if new chain is concatenated to old model"
        chain_id_increment=1
        .type=int
        .help="Cain id increment"
      }
    }
  }
  output
  .help="Output files"
  {
    logfile=xmanip.log
    .type=str
    .help="Logfile"
    hklout=xmanip.mtz
    .type=str
    .help="Ouptut miller indices and data"
    xyzout=xmanip.pdb
    .type=str
    .help="output PDB file"
  }
}
""", process_includes=True)

def print_help(command_name):
  print("""
#phil __OFF__
\t\t%s

A program for the manipulation of xray data objects (coordinates and reflection files).

The keywords are sumarized below:

#phil __ON__
xmanip {
  input {
    unit_cell = None
    space_group = None
    xray_data {
      file_name = None
      labels = None
      label_appendix = None
      name = None
      write_out = None
    }
    model {
      file_name = None
    }
  }
  parameters {
    action = reindex manipulate_pdb *manipulate_miller
    reindex {
      standard_laws = niggli *reference_setting invert user_supplied
      user_supplied_law = "h,k,l"
    }
    manipulate_miller {
      task = get_dano get_diso lsq_scale sfcalc *custom None
      output_label_root = "FMODEL"
      get_dano {
        input_data = None
      }
      get_diso {
        native = None
        derivative = None
        use_intensities = True
        use_weights = True
        scale_weight = True
      }
      lsq_scale {
        input_data_1 = None
        input_data_2 = None
        use_intensities = True
        use_weights = True
        scale_weight = True
      }
      sfcalc {
        fobs = None
        output = *2mFo-DFc mFo-DFc complex_fcalc abs_fcalc intensities
        use_bulk_and_scale = *as_estimated user_upplied
        bulk_and_scale_parameters {
          d_min = 2
          overall {
            b_cart {
              b_11 = 0
              b_22 = 0
              b_33 = 0
              b_12 = 0
              b_13 = 0
              b_23 = 0
            }
            k_overall = 0.1
          }
          solvent {
            k_sol = 0.3
            b_sol = 56
          }
        }
      }
      custom{
        code = print("hello world", file=out)
      }
    }
    manipulate_pdb{
      task = apply_operator *set_b
      apply_operator{
        operator = "x,y,z"
        invert=False
        concatenate_model=False
        chain_id_increment=1
      }
      set_b{
        b_iso = 30
      }
    }
  }
  output {
    logfile = "xmanip.log"
    hklout = "xmanip.mtz"
    xyzout = "xmanip.pdb"
  }
}
#phil __END__


Further details can be found in the documentation.

  """%(command_name))

def run(args, command_name="phenix.xmanip"):
  if (len(args) == 0 or "--help" in args or "--h" in args or "-h" in args):
    print_help(command_name=command_name)
  else:
    log = multi_out()
    if (not "--quiet" in args):
      log.register(label="stdout", file_object=sys.stdout)
    string_buffer = StringIO()
    string_buffer_plots = StringIO()
    log.register(label="log_buffer", file_object=string_buffer)

    phil_objects = []
    argument_interpreter = master_params.command_line_argument_interpreter(
      home_scope="map_coefs")

    print("#phil __OFF__", file=log)
    print("==========================", file=log)
    print("          XMANIP          ", file=log)
    print("reindexing and other tasks", file=log)
    print("==========================", file=log)
    print(file=log)


    for arg in args:
      command_line_params = None
      arg_is_processed = False
      # is it a file?
      if (os.path.isfile(arg)): ## is this a file name?
        # check if it is a phil file
        try:
          command_line_params = iotbx.phil.parse(file_name=arg)
          if command_line_params is not None:
            phil_objects.append(command_line_params)
            arg_is_processed = True
        except KeyboardInterrupt: raise
        except Exception : pass
      else:
        try:
          command_line_params = argument_interpreter.process(arg=arg)
          if command_line_params is not None:
            phil_objects.append(command_line_params)
            arg_is_processed = True
        except KeyboardInterrupt: raise
        except Exception : pass

      if not arg_is_processed:
        print("##----------------------------------------------##", file=log)
        print("## Unknown file or keyword:", arg, file=log)
        print("##----------------------------------------------##", file=log)
        print(file=log)
        raise Sorry("Unknown file or keyword: %s" % arg)

    effective_params = master_params.fetch(sources=phil_objects)
    params = effective_params.extract()

    # now get the unit cell from the files
    hkl_xs = []
    pdb_xs = None

    #multiple file names are allowed
    for xray_data in params.xmanip.input.xray_data:
      if xray_data.file_name is not None:
        hkl_xs.append( crystal_symmetry_from_any.extract_from(
           file_name=xray_data.file_name) )

    if params.xmanip.input.model.file_name is not None:
      pdb_xs = crystal_symmetry_from_any.extract_from(
        file_name=params.xmanip.input.model.file_name)

    phil_xs = crystal.symmetry(
      unit_cell=params.xmanip.input.unit_cell,
      space_group_info=params.xmanip.input.space_group  )

    combined_xs = crystal.select_crystal_symmetry(
      None,phil_xs, [pdb_xs],hkl_xs)
    if combined_xs is not None:
      # inject the unit cell and symmetry in the phil scope please
      params.xmanip.input.unit_cell = combined_xs.unit_cell()
      params.xmanip.input.space_group = \
        sgtbx.space_group_info( group = combined_xs.space_group() )

    print("#phil __ON__", file=log)
    new_params =  master_params.format(python_object=params)
    new_params.show(out=log)
    print("#phil __END__", file=log)

    if params.xmanip.input.unit_cell is None:
      raise Sorry("unit cell not specified")
    if params.xmanip.input.space_group is None:
      raise Sorry("space group not specified")

    #-----------------------------------------------------------
    #
    # step 1: read in the reflection file
    #

    miller_arrays = []
    labels = []
    label_appendix = []
    write_it = []
    names = {}

    if len(params.xmanip.input.xray_data)>0:

      phil_xs = crystal.symmetry(
        unit_cell=params.xmanip.input.unit_cell,
        space_group_info=params.xmanip.input.space_group  )

      xray_data_server =  reflection_file_utils.reflection_file_server(
        crystal_symmetry = phil_xs,
        force_symmetry = True,
        reflection_files=[])

      count=0
      for xray_data in params.xmanip.input.xray_data:
        if xray_data.file_name is not None:
          miller_array = None
          miller_array = read_data(xray_data.file_name,
                                   xray_data.labels,
                                   phil_xs)
          print(file=log)
          print("Summary info of observed data", file=log)
          print("=============================", file=log)
          if miller_array is None:
            raise Sorry("Failed to read data. see errors above" )
          miller_array.show_summary(f=log)
          print(file=log)

          miller_arrays.append( miller_array )
          labels.append( miller_array.info().labels )
          label_appendix.append( xray_data.label_appendix )

          this_name = "COL_"+str(count)
          if xray_data.name is not None:
            this_name = xray_data.name
          #check if this name is allready used
          if this_name in names:
            raise Sorry( "Non unique dataset name. Please change the input script" )
          names.update( {this_name:count} )
          count += 1

          write_it.append( xray_data.write_out)

      output_label_root = construct_output_labels( labels, label_appendix )
      for ii in range(len(labels)):
        test=0
        for jj in range( ii+1,len(labels) ):
          for lab_name1, lab_name2 in zip(labels[ii],labels[jj]):
            if lab_name1==lab_name2:
              test+=1
          if test == 2:
            print("\n***** You are trying to import the data with label(s) %s more then one time. ***** \n"%(str(labels[ii])), file=log)
      for ii in range(len(output_label_root)):
        for jj in range(ii+1,len(output_label_root)):
          if output_label_root[ii]==output_label_root[jj]:
            if write_it[ii]:
              if write_it[jj]:
                print("Output label roots:", file=log)
                print(output_label_root, file=log)
                raise Sorry( "Output labels are not unique. Modify input." )



    #----------------------------------------------------------------
    # Step 2: get an xray structure from the PDB file
    #
    pdb_model = None
    model = None
    if params.xmanip.input.model.file_name is not None:
      pdb_model = iotbx.pdb.input(
        file_name=params.xmanip.input.model.file_name)
      model = pdb_model.xray_structure_simple(crystal_symmetry=phil_xs)
      print("Atomic model summary", file=log)
      print("====================", file=log)
      model.show_summary(f=log)
      print(file=log)


    write_miller_array = False
    write_pdb_file = False
    # define some output holder thingamebobs
    new_miller_arrays = []
    new_model = None

    #manipulate miller arrays
    if params.xmanip.parameters.action == "manipulate_miller":
      write_miller_array = True
      new_miller = xmanip_tasks.manipulate_miller(names,
                                                  miller_arrays,
                                                  model,
                                                  params.xmanip.parameters.manipulate_miller,
                                                  log )
      miller_arrays.append( new_miller )
      # not very smart to rely here on a phil defintion defined in another file
      tmp_root = params.xmanip.parameters.manipulate_miller.output_label_root
      if tmp_root is None:
        tmp_root = "UNSPECIFIED"
      output_label_root.append( tmp_root )
      write_it.append(True)




    if params.xmanip.parameters.action=="reindex":
      write_miller_array = True
      #----------------------------------------------------------------
      # step 3: get the reindex laws
      phil_xs.show_summary()
      to_niggli    = phil_xs.change_of_basis_op_to_niggli_cell()
      to_reference = phil_xs.change_of_basis_op_to_reference_setting()
      to_inverse   = phil_xs.change_of_basis_op_to_inverse_hand()
      to_primitive = phil_xs.change_of_basis_op_to_primitive_setting()
      cb_op = None
      if (params.xmanip.parameters.reindex.standard_laws == "niggli"):
        cb_op = to_niggli
      if (params.xmanip.parameters.reindex.standard_laws == "reference_setting"):
        cb_op = to_reference
      if (params.xmanip.parameters.reindex.standard_laws == "invert"):
        cb_op = to_inverse
      if (params.xmanip.parameters.reindex.standard_laws == "user_supplied"):
        cb_op = sgtbx.change_of_basis_op( params.xmanip.parameters.reindex.user_supplied_law )
      if (params.xmanip.parameters.reindex.standard_laws == "primitive_setting"):
        cb_op = to_primitive


      if cb_op is None:
        raise Sorry("No change of basis operation is supplied.")

      print("Supplied reindexing law:", file=log)
      print("========================", file=log)
      print("hkl notation: ", cb_op.as_hkl(), file=log)
      print("xyz notation: ", cb_op.as_xyz(), file=log)
      print("abc notation: ", cb_op.as_abc(), file=log)
      #----------------------------------------------------------------
      # step 4: do the reindexing
      #
      # step 4a: first do the miller array object
      #new_miller_arrays = []
      for miller_array in miller_arrays:
        new_miller_array = None
        if miller_array is not None:
          new_miller_array = miller_array.change_basis( cb_op )
          new_miller_arrays.append( new_miller_array )
      #
      # step 4b: the xray structure
      if pdb_model is not None:
        write_pdb_file=True
        new_model = model.change_basis( cb_op )


    if write_miller_array:
      if len(new_miller_arrays)==0:
        new_miller_arrays = miller_arrays
      #----------------------------------------------------------------
      print(file=log)
      print("The data has been reindexed/manipulated", file=log)
      print("--------------------------------------", file=log)
      print(file=log)
      print("Writing output files....", file=log)

      mtz_dataset=None
      if len(new_miller_arrays)>0:
        first=0
        for item in range(len(write_it)):
          if write_it[item]:
            first=item
            if new_miller_arrays[ first ] is not None:
              break

        if new_miller_arrays[first] is not None:
          tmp = new_miller_arrays[first].map_to_asu()
          mtz_dataset = tmp.as_mtz_dataset(
            column_root_label=output_label_root[first])

      if mtz_dataset is not None:
        for miller_array, new_root in zip(new_miller_arrays[first+1:],
                                          output_label_root[first+1:]):
          if miller_array is not None:
            mtz_dataset = mtz_dataset.add_miller_array(
              miller_array = miller_array,
              column_root_label = new_root)

        print("Writing mtz file with name %s"%(params.xmanip.output.hklout), file=log)
        mtz_dataset.mtz_object().write(
          file_name=params.xmanip.output.hklout)

      #step 5b: write the new pdb file
      if new_model is not None:
        pdb_file = open( params.xmanip.output.xyzout, 'w')
        print("Wring pdb file to: %s"%(params.xmanip.output.xyzout), file=log)
        write_as_pdb_file(
          input_pdb = pdb_model,
          input_xray_structure = new_model,
          out = pdb_file,
          chain_id_increment= 0,
          additional_remark = "Generated by %s" % command_name)

        pdb_file.close()
      if ( [miller_array,new_model]).count(None)==2:
        print("No input reflection of coordinate files have been given", file=log)

    if params.xmanip.parameters.action=="manipulate_pdb":
      if params.xmanip.parameters.manipulate_pdb.task == "apply_operator":
        rt_mx = None
        if params.xmanip.parameters.manipulate_pdb.apply_operator.standard_operators == "user_supplied_operator":
          rt_mx = sgtbx.rt_mx(
            params.xmanip.parameters.manipulate_pdb.apply_operator.user_supplied_operator,t_den=12*8 )
          print("Applied operator : ", rt_mx.as_xyz(), file=log)
        if params.xmanip.parameters.manipulate_pdb.apply_operator.standard_operators == \
             "user_supplied_cartesian_rotation_matrix":
          rt = params.xmanip.parameters.manipulate_pdb.apply_operator.user_supplied_cartesian_rotation_matrix
          tmp_r=None
          tmp_t=None
          if "," in rt.r:
            tmp_r = rt.r.split(',')
          else:
            tmp_r = rt.r.split(' ')
          if "," in rt.r:
            tmp_t = rt.t.split(',')
          else:
            tmp_t = rt.t.split(' ')
          tmp_tmp_r=[]
          tmp_tmp_t=[]
          for item in tmp_r:
            tmp_tmp_r.append( float(item) )
          if len(tmp_tmp_r)!=9:
            raise Sorry("Invalid rotation matrix. Please check input: %s"%(rt.r) )
          for item in tmp_t:
            tmp_tmp_t.append( float(item) )
          if len(tmp_tmp_t)!=3:
            raise Sorry("Invalid translational vector. Please check input: %s"%(rt.t) )
          tmp_tmp_t = (tmp_tmp_t)
          rt_mx = quick_rt_mx(tmp_tmp_r, tmp_tmp_t)
          print("User supplied cartesian matrix and vector: ", file=log)
          rt_mx.show()
          o = matrix.sqr(model.unit_cell().orthogonalization_matrix())
          tmp_r = o.inverse()*rt_mx.r()*o
          tmp_t = o.inverse()*matrix.col(list(rt_mx.t()))
          print(file=log)
          print("Operator in fractional coordinates: ", file=log)
          rt_mx = quick_rt_mx(r=tmp_r.as_float(), t=list(tmp_t))
          rt_mx.show(out=log)
          print(file=log)


        if params.xmanip.parameters.manipulate_pdb.apply_operator.invert:
          rt_mx = rt_mx.inverse()
          print(file=log)
          print("Taking inverse of given operator", file=log)
          print(file=log)

        sites = model.sites_frac()
        new_sites = flex.vec3_double()
        for site in sites:
          new_site = rt_mx.r()*matrix.col(site)
          new_site = flex.double(new_site)+flex.double( rt_mx.t().as_double() )
          new_sites.append( tuple(new_site) )
        new_model = model.deep_copy_scatterers()

        new_model.set_sites_frac( new_sites )
        # write the new [pdb file please
        pdb_file = open( params.xmanip.output.xyzout, 'w')
        print("Wring pdb file to: %s"%(params.xmanip.output.xyzout), file=log)
        if params.xmanip.parameters.manipulate_pdb.apply_operator.concatenate_model:
          write_as_pdb_file( input_pdb = pdb_model,
                             input_xray_structure = model,
                             out = pdb_file,
                             chain_id_increment = 0,
                             additional_remark = None,
                             print_cryst_and_scale=True )

        write_as_pdb_file( input_pdb = pdb_model,
                           input_xray_structure = new_model,
                           out = pdb_file,
                           chain_id_increment = params.xmanip.parameters.manipulate_pdb.apply_operator.chain_id_increment,
                           additional_remark = None,
                           print_cryst_and_scale=False )

        pdb_file.close()

      if params.xmanip.parameters.manipulate_pdb.task =="set_b":
        #rest all the b values
        if params.xmanip.parameters.manipulate_pdb.set_b:
          b_iso = params.xmanip.parameters.manipulate_pdb.b_iso
          new_model = model.set_b_iso( value = b_iso )
          print(file=log)
          print("All B-values have been set to %5.3f"%(b_iso), file=log)
          print("Writing PDB file %s"%(params.xmanip.output.xyzout), file=log)
          print(file=log)

        pdb_file = open( params.xmanip.output.xyzout, 'w')
        write_as_pdb_file( input_pdb = pdb_model,
                           input_xray_structure = new_model,
                           out = pdb_file,
                           chain_id_increment = 0,
                           additional_remark = None,
                           print_cryst_and_scale=True)
        pdb_file.close()

    #write the logfile
    logger = open( params.xmanip.output.logfile, 'w')
    print("Writing log file with name %s  "%(params.xmanip.output.logfile), file=log)
    print(string_buffer.getvalue()[0:len(string_buffer.getvalue())-1], file=logger) #avoid a newline at the end ...
    logger.close()
