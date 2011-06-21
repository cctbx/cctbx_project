## Copyright 2002 by PyMMLib Development Group (see AUTHORS file)
## This code is part of the PyMMLib distrobution and governed by
## its license.  Please see the LICENSE file that should have been
## included as part of this package.
#
# http://cvs.sourceforge.net/viewcvs.py/pymmlib/pymmlib/LICENSE.txt?view=markup
#
# based on PyMMLib 0.5
# minor changes by R.W. Grosse-Kunstleve (communicated to PyMMLib authors)
#
"""mmCIF file and mmCIF dictionary parser.  Files are parsed
into a set of data structures where they can be further processed.  The
data structures can also be constructed and written back out as mmCIF.
A CIF dictionary parser is also included as a specialized version of the
mmCIF parser.
"""
import string
import re
import copy
import sys
import os
from types import IntType, StringType

def OpenFile(path, mode):
    """Right now this only supports opening GZip'ed files, in the future
    it might be extended for URLs.
    """
    ## if path is not a string, assume it is a file object and
    ## return it
    if type(path) != StringType:
        return path
    (base, ext) = os.path.splitext(path)
    if ext == ".gz":
        return gzip.open(path, mode)
    return open(path, mode)

## errors
mmCIFError = "mmCIFError"


##
## DATA STRUCTURES FOR HOLDING CIF INFORMATION
##
## mmCIF files are parsed into:
##         mmCIFFile -> [mmCIFData] -> [mmCIFTable] -> [mmCIFRow]
##
## mmCIF dictionaries are parsed into:
##         mmCIFDictionary -> [mmCIFData] -> [mmCIFTable] -> [mmCIFRow]
##

mmCIFError = "mmCIFError"
MAX_LINE = 80

class mmCIFRow(dict):
    """Contains one row of data.  In a mmCIF file, this is one complete
    set of data found under a section.  The data can be accessed by using
    the column names as class attributes.
    """
    __slots__ = ["table"]

    def __eq__(self, other):
        return id(self) == id(other)

    def __deepcopy__(self, memo):
        return mmCIFRow(self)

    def mget(self, *keys):
        """Return the fist value found for the given keys in the argument
        list.
        """
        for key in keys:
            try:
                return self[key]
            except KeyError:
                continue


class mmCIFTable(list):
    """Contains columns and rows of data for a mmCIF section.  Rows of data
    are stored as mmCIFRow classes.
    """
    __slots__ = ["name", "columns", "data"]

    def __init__(self, name, columns = None):
        list.__init__(self)
        self.name = name
        self.columns = columns or []

    def __deepcopy__(self, memo):
        table = mmCIFTable(self.name, self.columns[:])
        for row in self:
            table.append(copy.deepcopy(row, memo))
        return table

    def __eq__(self, other):
        return id(self) == id(other)

    def __getitem__(self, x):
        """Retrieves mmCIFRow at index x from the table if the argument is
        a integer.  If the argument is a string, then the data from the
        first row is returned.
        """
        if type(x) == IntType:
            return list.__getitem__(self, x)

        elif type(x) == StringType:
            try:
                return self[0][x]
            except IndexError:
                raise KeyError
            except KeyError:
                raise KeyError

        raise TypeError, x

    def __setitem__(self, i, row):
        assert isinstance(row, mmCIFRow)
        row.table = self
        list.__setitem__(self, i, row)

    def __delitem__(self, i):
        assert isinstance(row, mmCIFRow)
        self.remove(self[i])

    def append(self, row):
        assert isinstance(row, mmCIFRow)
        row.table = self
        list.append(self, row)

    def insert(self, i, row):
        assert isinstance(row, mmCIFRow)
        row.table = self
        list.insert(self, i, row)

    def remove(self, row):
        assert isinstance(row, mmCIFRow)
        del row.table
        list.remove(self, row)

    def autoset_columns(self):
        """Iterates through all rows in self, and forms a list of all
        unique column names, then sets the self.columns to that list.
        """
        self.columns = []
        for row in self:
            for col in row.keys():
                if col not in self.columns:
                    self.columns.append(col)
        self.columns.sort()

    def get_row(self, *args):
        """Preforms a SQL-like 'AND' select aginst all the rows in the table,
        and returns the first matching row found.  The arguments are a
        variable list of tuples of the form:
          (<column-name>, <column-value>)
        For example:
          ger_row(('atom_id','CA'),('entity_id', '1'))
        returns the first matching row with atom_id==1 and entity_id==1.
        """
        def cmp(r):
            try:
                for (k, v) in args:
                    if r[k] != v:
                        return False
            except KeyError:
                return False
            return True

        for row in self:
            if cmp(row):
                return row
        return None

    def iter_rows(self, *args):
        """This is the same as get_row, but it iterates over all matching
        rows in the table.
        """
        def cmp(r):
            try:
                for (k, v) in args:
                    if r[k] != v:
                        return False
            except KeyError:
                return False
            return True

        for row in self:
            if cmp(row):
                yield row

    def row_index_dict(self, key):
        """Return a dictionary mapping the value of the row's value in
        column 'key' to the row itself.  If there are multiple rows with
        the same key value, they will be overwritten with the last found
        row.
        """
        dictx = {}
        for row in self:
            try:
                dictx[row[key]] = row
            except KeyError:
                pass
        return dictx

    def debug(self):
        print "mmCIFTable::%s" % (self.name)
        for row in self:
            for col in self.columns:
                print "%s=%s" % (col, row.get(col))[:80]
            print "---"


class mmCIFData(list):
    """Contains all information found under a data_ block in a mmCIF file.
    mmCIF files are represented differently here than their file format
    would suggest.  Since a mmCIF file is more-or-less a SQL database dump,
    the files are represented here with their sections as "Tables" and
    their subsections as "Columns".  The data is stored in "Rows".
    """
    __slots__ = ["name", "file"]

    def __init__(self, name):
        list.__init__(self)
        self.name = name

    def __deepcopy__(self, memo):
        data = mmCIFData(self.name)
        for table in self:
            data.append(copy.deepcopy(table, memo))
        return data

    def __eq__(self, other):
        return id(self) == id(other)

    def __getitem__(self, x):
        if type(x) == IntType:
            return list.__getitem__(self, x)

        elif type(x) == StringType:
            for ctable in self:
                if x == ctable.name:
                    return ctable
            raise KeyError, x

        raise TypeError, x

    def __delitem__(self, x):
        """Remove a mmCIFTable by index or table name.
        """
        self.remove(self[x])

    def append(self, table):
        """Append a mmCIFTable.  This will trigger the removal of any
        table with the same name.
        """
        assert isinstance(table, mmCIFTable)
        try:
            del self[table.name]
        except KeyError:
            pass
        table.data = self
        list.append(self, table)

    def insert(self, i, table):
        assert isinstance(table, mmCIFTable)
        try:
            del self[table.name]
        except KeyError:
            pass
        table.data = self
        list.insert(self, i, table)

    def remove(self, table):
        assert isinstance(table, mmCIFTable)
        del table.data
        list.remove(self, table)

    def has_key(self, x):
        try:
            self[x]
        except KeyError:
            return False
        else:
            return True

    def get(self, x, default = None):
        try:
            return self[x]
        except KeyError:
            return default

    def get_table(self, name):
        """Looks up and returns a stored mmCIFTable class by it's name.  This
        name is the section key in the mmCIF file.
        """
        try:
            return self[name]
        except KeyError:
            return None
        except IndexError:
            return None

    def debug(self):
        print "mmCIFData::%s" % (self._name)
        for ctable in self:
            ctable.debug()


class mmCIFSave(mmCIFData):
    """Class to store data from mmCIF dictionary save_ blocks.  I treat
    them as non-nested sections along with data_ sections.  This may
    not be correct.
    """
    pass


class mmCIFFile(list):
    """Class representing a mmCIF files.
    """
    def __deepcopy__(self, memo):
        file = mmCIFFile()
        for data in self:
            file.append(copy.deepcopy(data, memo))
        return file

    def __eq__(self, other):
        return id(self) == id(other)

    def __getitem__(self, x):
        """Retrieve a mmCIFData object by index or name.
        """
        if type(x) == IntType:
            return list.__getitem__(self, x)

        elif type(x) == StringType:
            for cdata in self:
                if cdata.name == x:
                    return cdata
            raise KeyError, x

        raise TypeError, x

    def __delitem__(self, x):
        """Remove a mmCIFData by index or data name.  Raises IndexError
        or KeyError if the mmCIFData object is not found, the error raised
        depends on the argument type.
        """
        self.remove(self[x])

    def append(self, cdata):
        """Append a mmCIFData object.  This will trigger the removal of any
        mmCIFData object in the file with the same name.
        """
        assert isinstance(cdata, mmCIFData)
        try:
            del self[cdata.name]
        except KeyError:
            pass
        cdata.file = self
        list.append(self, cdata)

    def insert(self, i, cdata):
        assert isinstance(cdata, mmCIFData)
        try:
            del self[cdata.name]
        except KeyError:
            pass
        cdata.file = self
        list.insert(self, i, cdata)

    def has_key(self, x):
        for cdata in self:
            if cdata.name == x:
                return True
        return False

    def get(self, x, default = None):
        try:
            return self[x]
        except KeyError:
            return default

    def load_file(self, fil, update_cb = None, strict = True):
        """Load and append the mmCIF data from file object fil into self.
        """
        if (isinstance(fil, str)):
          fil = OpenFile(fil, "r")
        mmCIFFileParser(strict=strict).parse_file(fil, self, update_cb)

    def save_file(self, fil):
        fil = OpenFile(fil, "w")
        mmCIFFileWriter().write_file(fil, self)

    def get_data(self, name):
        try:
            return self[name]
        except KeyError:
            return None
        except IndexError:
            return None

    def debug(self):
        print "mmCIFFile"
        for cdata in self:
            cdata.debug()


class mmCIFDictionary(mmCIFFile):
    """Class representing a mmCIF dictionary.  The constructor of this class
    takes two arguments.  The first is the string path for the file, or
    alternativly a file object.
    """
    pass


##
## FILE PARSERS/WRITERS
##


class mmCIFFileParser(object):
    """Stateful parser which uses the mmCIFElementFile tokenizer to read
    a mmCIF file and convert it into the mmCIFData/mmCIFTable/mmCIFRow
    data hierarchy.
    """
    def __init__(self, strict = True):
        self._strict = strict

    def parse_file(self, fil, cif_file, update_cb = None):
        self.update_cb = update_cb
        self.line_number = 0

        token_iter = self.gen_token_iter(fil)

        try:
            self.parse(token_iter, cif_file)
        except StopIteration:
            pass
        else:
            raise mmCIFError

    def syntax_error(self, err):
        err = "[line %d] %s" % (self.line_number, err)
        raise mmCIFError, err

    def parse(self, token_iter, cif_file):
        cif_table_cache = {}
        cif_data        = None
        cif_table       = None
        cif_row         = None
        state           = ""

        tblx,colx,strx,tokx = token_iter.next()

        while 1:
            if tblx != None:
                state = "RD_SINGLE"
            elif tokx != None:
                if tokx == "loop_":
                    state = "RD_LOOP"
                elif tokx.startswith("data_"):
                    state = "RD_DATA"
                elif tokx.startswith("save_"):
                    state = "RD_SAVE"
                elif (self._strict):
                    self.syntax_error("bad token #1: "+str(tokx))
                else:
                    ## ignore input until a data_ or save_ token is found
                    while 1:
                        tblx,colx,strx,tokx = token_iter.next()
                        if (tokx is None): continue
                        if tokx.startswith("data_"):
                            state = "RD_DATA"
                            break
                        elif tokx.startswith("save_"):
                            state = "RD_SAVE"
                            break
            else:
                self.syntax_error("bad token #2")
                return


            if state == "RD_SINGLE":
                try:
                    cif_table = cif_table_cache[tblx]
                except KeyError:
                    cif_table = cif_table_cache[tblx] = mmCIFTable(tblx)

                    try:
                        cif_data.append(cif_table)
                    except AttributeError:
                        self.syntax_error(
                            "section not contained in data_ block")
                        return

                    cif_row = mmCIFRow()
                    cif_table.append(cif_row)
                else:
                    try:
                        cif_row = cif_table[0]
                    except IndexError:
                        self.syntax_error("bad token #3")
                        return

                ## check for duplicate entries
                if colx in cif_table.columns:
                    self.syntax_error("redefined subsection (column)")
                    return
                else:
                    cif_table.columns.append(colx)

                x,x,strx,tokx = token_iter.next()

                if tokx != None:
                    if tokx == ".":
                        cif_row[colx] = ""
                    else:
                        cif_row[colx] = tokx
                elif strx != None:
                    cif_row[colx] = strx
                else:
                    self.syntax_error("bad token #4")

                tblx,colx,strx,tokx = token_iter.next()
                continue

            elif state == "RD_LOOP":
                tblx,colx,strx,tokx = token_iter.next()

                if tblx == None or colx == None:
                    self.syntax_error("bad token #5")
                    return

                if cif_table_cache.has_key(tblx):
                    self.syntax_error("_loop section duplication")
                    return

                cif_table = mmCIFTable(tblx)

                try:
                    cif_data.append(cif_table)
                except AttributeError:
                    self.syntax_error(
                        "_loop section not contained in data_ block")
                    return

                cif_table.columns.append(colx)

                while 1:
                    tblx,colx,strx,tokx = token_iter.next()
                    if tblx == None:
                        break
                    if tblx != cif_table.name:
                        self.syntax_error("changed section names in _loop")
                        return
                    cif_table.columns.append(colx)

                if (not self._strict):
                    if tokx in ["loop_", "data_", "save_"]:
                        continue

                while 1:
                    cif_row = mmCIFRow()
                    cif_table.append(cif_row)

                    for col in cif_table.columns:
                        if tokx in ["loop_", "data_", "save_"]:
                          self.syntax_error(
                            "unexpected reserved word in loop: %s" % tokx)
                        if tokx != None:
                            if tokx == ".":
                                cif_row[col] = ""
                            else:
                                cif_row[col] = tokx
                        elif strx != None:
                            cif_row[col] = strx

                        tblx,colx,strx,tokx = token_iter.next()

                    ## the loop ends when one of these conditions is met
                    if tblx != None:
                        break
                    if tokx == None:
                        continue
                    if tokx == "loop_":
                        break
                    if tokx.startswith("data_"):
                        break
                    if tokx.startswith("save_"):
                        break

                continue

            elif state == "RD_DATA":
                cif_data = mmCIFData(tokx[5:])
                cif_file.append(cif_data)
                cif_table_cache = {}
                cif_table = None

                tblx,colx,strx,tokx = token_iter.next()

            elif state == "RD_SAVE":
                cif_data = mmCIFSave(tokx[5:])
                cif_file.append(cif_data)
                cif_table_cache = {}
                cif_table = None

                tblx,colx,strx,tokx = token_iter.next()


    def gen_token_iter(self, fil):
        re_tok = re.compile(
            r"(?:"

             "(?:_(.+?)[.](\S+))"         "|"  # _section.subsection

             "(?:['](.*?)(?:[']\s|[']$))" "|"  # quoted strings

             "(?:\s#.*$)"                 "|"  # comments

             "(\S+)"                           # unquoted tokens

             ")")

        ## get file size for update callbacks
        percent_done = 0
        fil_read_bytes = 0

        ## some file objects do not support seek/tell
        if hasattr(fil, "seek") and hasattr(fil, "tell"):
            try:
                fil.seek(0, 2)
                fil_size_bytes = fil.tell()
                fil.seek(0, 0)
            except:
                # this is a adverage file size ;)
                fil_size_bytes = 1304189
        else:
            fil_size_bytes = 1304189

        file_iter = iter(fil)
        while 1:
            try:
                ln = file_iter.next()
            except StopIteration:
                break
            else:
                self.line_number += 1
                fil_read_bytes   += len(ln)

            ## call update callback
            if self.update_cb != None:
                pdone = (fil_read_bytes * 100)/fil_size_bytes
                if pdone != percent_done and pdone <= 100:
                    percent_done = pdone
                    self.update_cb(percent_done)

            ## skip comments
            if ln.startswith("#"):
                continue

            ## semi-colen multi-line strings
            if ln.startswith(";"):
                x = ln[1:]
                while 1:
                    try:
                        ln = file_iter.next()
                    except StopIteration:
                        break
                    else:
                        self.line_number += 1

                    if ln.startswith(";"):
                        break
                    x += ln

                x = x.rstrip()
                yield (None, None, x, None)
                continue

            ## split line into tokens
            tok_iter = re_tok.finditer(ln)
            for tokm in tok_iter:
                if (tokm.groups() != (None,None,None,None)):
                    yield tokm.groups()


class mmCIFFileWriter(object):
    """Writes out a mmCIF file using the data in the mmCIFData list.
    """
    def write_file(self, fil, cif_data_list):
        self.fil = fil

        ## constant controlls the spacing between columns
        self.SPACING = 2

        ## iterate through the data sections and write them
        ## out to the file
        for cif_data in cif_data_list:
            self.cif_data = cif_data
            self.write_cif_data()

    def write(self, x):
        self.fil.write(x)

    def writeln(self, x = ""):
        self.fil.write(x+"\n")

    def write_mstring(self, mstring):
        self.write(self.form_mstring(mstring))

    def form_mstring(self, mstring):
        strx = ";"

        lw = MAX_LINE - 2

        for x in mstring.split("\n"):
            if x == "":
                strx += "\n"
                continue

            while len(x) > 0:
                x1 = x[:lw]
                x  = x[lw:]
                strx += x1 + "\n"

        strx += ";\n"
        return strx

    def fix(self, x):
        if type(x) != StringType:
            return str(x)
        if x == "":
            return "."
        return x

    def data_type(self, x):
        """Analyze x and return its type: token, qstring, mstring
        """
        if type(x) != StringType:
            x = str(x)
            return x, "token"

        if x == "":
            return ".", "token"

        if x.find("\n") != -1:
            return x, "mstring"

        if x.find(" ") != -1 or x.find("\t") != -1:
            if len(x) > MAX_LINE-2:
                return x, "mstring"
            if x.find("' ") != -1:
                return x, "mstring"
            return x, "qstring"

        if len(x) < MAX_LINE:
            return x, "token"
        else:
            return x, "mstring"

    def write_cif_data(self):
        if isinstance(self.cif_data, mmCIFSave):
            self.writeln("save_%s" % self.cif_data.name)
        else:
            self.writeln("data_%s" % self.cif_data.name)
        self.writeln("#")

        for cif_table in self.cif_data:
            ## ignore tables without data rows
            if len(cif_table) == 0:
                continue

            ## special handling for tables with one row of data
            elif len(cif_table) == 1:
                self.write_one_row_table(cif_table)

            ## _loop tables
            elif len(cif_table) > 1 and len(cif_table.columns) > 0:
                self.write_multi_row_table(cif_table)

            else:
                print "wtf?",cif_table
                sys.exit(1)

            self.writeln("#")

    def write_one_row_table(self, cif_table):
        row = cif_table[0]

        ## determine max key length for formatting output
        kmax  = 0
        table_len = len(cif_table.name) + 2
        for col in cif_table.columns:
            klen = table_len + len(col)
            assert klen < MAX_LINE
            kmax = max(kmax, klen)

        ## we need a space after the tag
        kmax += self.SPACING
        vmax  = MAX_LINE - kmax - 1

        ## write out the keys and values
        strx = ""

        for col in cif_table.columns:
            strx = "_%s.%s" % (cif_table.name, col)
            strx = strx.ljust(kmax)

            try:
                x0 = row[col]
            except KeyError:
                x = "?"
                dtype = "token"
            else:
                x, dtype = self.data_type(x0)

            if dtype == "token":
                if len(x) > vmax:
                    strx += "\n"
                strx += "%s\n" % (x)
                self.write(strx)

            elif dtype == "qstring":
                if len(x) > vmax:
                    strx += "\n"
                strx += "'%s'\n" % (x)
                self.write(strx)

            elif dtype == "mstring":
                strx += "\n"
                self.write(strx)
                self.write_mstring(x)

    def write_multi_row_table(self, cif_table):
        ## write the key description for the _loop
        self.writeln("loop_")
        for col in cif_table.columns:
            key = "_%s.%s" % (cif_table.name, col)
            assert len(key) < MAX_LINE
            self.writeln(key)

        col_len_map   = {}
        col_dtype_map = {}

        for row in cif_table:
            for col in cif_table.columns:
                ## get data and data type
                try:
                    x0 = row[col]
                except KeyError:
                    lenx  = 1
                    dtype = "token"
                else:
                    x, dtype = self.data_type(x0)

                    ## determine write length of data
                    if dtype == "token":
                        lenx = len(x)
                    elif dtype == "qstring":
                        lenx = len(x) + 2
                    else:
                        lenx = 0

                try:
                    col_dtype = col_dtype_map[col]
                except KeyError:
                    col_dtype_map[col] = dtype
                    col_len_map[col] = lenx
                    continue

                ## modify column data type if necessary
                if col_dtype != dtype:
                    if col_dtype == "mstring":
                        continue
                    elif (col_dtype == "qstring" or col_dtype == "token") and \
                         dtype == "mstring":
                        col_dtype_map[col] = "mstring"
                        continue
                    elif col_dtype == "token" and dtype == "qstring":
                        col_dtype_map[col] = "qstring"

                ## update the column charactor width if necessary
                if col_len_map[col] < lenx:
                    col_len_map[col] = lenx

        ## form a write list of the column names with values of None to
        ## indicate a newline
        wlist = []
        llen = 0
        for col in cif_table.columns:
            dtype = col_dtype_map[col]

            if dtype == "mstring":
                llen = 0
                wlist.append((None, None, None))
                wlist.append((col, dtype, None))
                continue

            lenx  = col_len_map[col]
            if llen == 0:
                llen = lenx
            else:
                llen += self.SPACING + lenx

            if llen > MAX_LINE-1:
                wlist.append((None, None, None))
                llen = lenx

            wlist.append((col, dtype, lenx))

        ## write out the data
        spacing = " " * self.SPACING
        add_space = False
        listx = []

        for row in cif_table:
            for (col, dtype, lenx) in wlist:

                if col == None:
                    add_space = False
                    listx.append("\n")
                    continue

                elif add_space == True:
                    add_space = False
                    listx.append(spacing)


                if dtype == "token":
                    try:
                        x = row[col]
                    except KeyError:
                        x = "?"
                    else:
                        if type(x) != StringType:
                            x = str(x)
                        if x == "":
                            x = "."

                    x = x.ljust(lenx)
                    listx.append(x)
                    add_space = True


                elif dtype == "qstring":
                    try:
                        x = "'%s'" % (row[col])
                    except KeyError:
                        x = "?".ljust(lenx)

                    x = x.ljust(lenx)
                    listx.append(x)
                    add_space = True


                elif dtype == "mstring":
                    try:
                        listx.append(self.form_mstring(row[col]))
                    except KeyError:
                        listx.append("?\n")
                    add_space = False


            add_space = False
            listx.append("\n")

            ## write out strx if it gets big to avoid using a lot of
            ## memory
            if len(listx) > 1024:
                self.write("".join(listx))
                listx = []

        ## write out the _loop section
        self.write("".join(listx))


mmCIFStandardColumnsMap = {
    "entry":        ["id"],

    "entity":       ["id", "type", "details"],

    "audit_author": ["name"],

    "cell":         ["entry_id", "length_a", "length_b", "length_c",
                     "angle_alpha", "angle_beta", "angle_gamma", "PDB_Z"],

    "symmetry":     ["entry_id", "space_group_name_H-M", "cell_setting",
                     "Int_Tables_number"],

    "atom_site":    ["group_PDB", "id", "type_symbol", "label_entity_id",
                     "Cartn_x", "Cartn_y", "Cartn_z",
                     "occupancy", "B_iso_or_equiv", "Cartn_x_esd",
                     "Cartn_y_esd", "Cartn_z_esd", "occupancy_esd",
                     "B_iso_or_equiv_esd",
                     "auth_asym_id",
                     "auth_seq_id",
                     "auth_comp_id",
                     "auth_alt_id",
                     "auth_atom_id"],

    "atom_site_anisotrop": [
                     "id", "type_symbol", "label_entity_id",
                     "U[1][1]", "U[1][2]", "U[1][3]", "U[2][2]",
                     "U[2][3]", "U[3][3]", "U[1][1]_esd", "U[1][2]_esd",
                     "U[1][3]_esd", "U[2][2]_esd", "U[2][3]_esd",
                     "U[3][3]_esd", "pdbx_auth_seq_id",
                     "pdbx_auth_comp_id", "pdbx_auth_asym_id",
                     "pdbx_auth_atom_id"]}


class mmCIFFileBuilder(object):
    """Builds a mmCIF file from a Structure object.
    """
    cifdb_omit_list = [
        "entity", "cell", "symmetry", "atom_site", "atom_site_anisotrop"]

    def __init__(self, struct, cif_file):
        self.struct = struct
        self.cif_data = mmCIFData("XXX")
        cif_file.append(self.cif_data)

        ## maps fragment -> entity_id
        self.entity_id_map = {}

        ## tables which are not generated from the structure hierarchy
        ## can be copied directly from the structure's cif database
        for table in self.struct.cifdb:
            if table.name not in self.cifdb_omit_list:
                new_table = copy.deepcopy(table)
                new_table.autoset_columns()
                self.cif_data.append(new_table)

        ## these tables need to be formed from the atom structure
        self.add__entry()
        self.add__entity()
        self.add__cell()
        self.add__symmetry()
        self.add__atom_site()

    def get_table(self, name, columns = None):
        try:
            return self.cif_data[name]
        except KeyError:
            pass

        if columns == None:
            columns = mmCIFStandardColumnsMap[name]

        table = mmCIFTable(name, columns[:])
        self.cif_data.append(table)
        return table

    def add__entry(self):
        """Add the _entry table.  If there is not entry ID, it defaults
        to XXX.
        """
        try:
            entry = self.cif_data["entry"]
        except KeyError:
            entry = self.get_table("entry")

        try:
            row0 = entry[0]
        except IndexError:
            row0 = mmCIFRow()
            entry.append(row0)

        if not row0.has_key("id"):
            row0["id"] = "XXX"

        self.cif_data.name = row0["id"]

    def add__entity(self):
        """Adds the entity table.  The entity names are faked here, since
        it's really not clear to me how the names are chosen by the PDB.
        """

        ## maps fragment -> entity_id
        entity = self.get_table("entity")

        ## I. entity.type == "polymer"

        ## first detect polymer chains
        ## map of entity::mmCIFRow -> sequence list
        es_list = []

        for chain in self.struct.iter_chains():

            ## if the chain is a bio-polymer, it is one entity; come up
            ## with a name from its sequence and add it to the
            ## entity map
            if not chain.has_standard_residues():
                continue

            sequence = chain.sequence or chain.calc_sequence()

            ## compare aginst previously calculated sequences to
            ## determine the correct entity_id
            entity_id = None

            for (row, seq) in es_list:
                if seq == sequence:
                    entity_id = row["id"]
                    break

            if entity_id == None:
                row = mmCIFRow()
                entity.append(row)

                row["id"] = entity.index(row) + 1
                row["type"] = "polymer"

                if self.struct.library.is_amino_acid(sequence[0]):
                    row["details"] = "%d residue polypeptide"%(len(sequence))
                elif self.struct.library.is_nucleic_acid(sequence[0]):
                    row["details"] = "%d residue DNA/RNA"%(len(sequence))

                entity_id = row["id"]
                es_list.append((row, sequence))

            for res in chain.iter_standard_residues():
                self.entity_id_map[res] = entity_id


        ## II. entity.type == "non-polymer" or "water"
        er_map = {}

        for chain in self.struct.iter_chains():
            for frag in chain.iter_non_standard_residues():

                ## already assigned a entity_id for this fragment_id
                if er_map.has_key(frag.res_name):
                    self.entity_id_map[frag] = er_map[frag.res_name]

                ## we need to assign a entity_id for this fragment_id
                ## and add a row for it in the entity table
                else:
                    row = mmCIFRow()
                    entity.append(row)

                    entity_id = row["id"] = entity.index(row) + 1

                    if frag.is_water():
                        row["type"] = "water"
                        row["details"] = ""
                    else:
                        row["type"] = "non-polymer"
                        row["details"] = frag.res_name

                    er_map[frag.res_name] = entity_id
                    self.entity_id_map[frag] = entity_id

    def add__cell(self):
        """Adds the _cell table.
        """
        if self.struct.unit_cell:
            unit_cell = self.struct.unit_cell
        else:
            return

        cell = self.get_table("cell")
        row = mmCIFRow()
        cell.append(row)

        row["entry_id"] = self.cif_data["entry"]["id"]
        row["length_a"] = unit_cell.a
        row["length_b"] = unit_cell.b
        row["length_c"] = unit_cell.c
        row["angle_alpha"] = unit_cell.calc_alpha_deg()
        row["angle_beta"] = unit_cell.calc_beta_deg()
        row["angle_gamma"] = unit_cell.calc_gamma_deg()

    def add__symmetry(self):
        """Adds the _symmetry table.
        """
        if self.struct.unit_cell and self.struct.unit_cell.space_group:
            space_group = self.struct.unit_cell.space_group
        else:
            return

        cell = self.get_table("symmetry")
        row = mmCIFRow()
        cell.append(row)

        row["entry_id"] = self.cif_data["entry"]["id"]
        row["space_group_name_H-M"] = space_group.pdb_name
        row["Int_Tables_number"] = space_group.number

    def add__atom_site(self):
        """Adds the _atom_site table.
        """
        atom_site = self.get_table("atom_site")

        orig_model = self.struct.default_model

        atom_id = 0

        for model_num in self.struct.model_list():
            self.struct.default_model = model_num

            for frag in self.struct.iter_fragments():
                for atm in frag.iter_all_alt_loc_atoms():

                    asrow = mmCIFRow()
                    atom_site.append(asrow)

                    atom_id += 1
                    asrow["id"] = atom_id

                    self.set_atom_site_row(asrow, atm)

        self.struct.default_model = orig_model

    def set_atom_site_row(self, asrow, atm):
        if atm.get_fragment().is_standard_residue():
            asrow["group_PDB"] = "ATOM"
        else:
            asrow["group_PDB"] = "HETATM"

        asrow["label_entity_id"] = self.entity_id_map[atm.get_fragment()]
        asrow["auth_atom_id"] = atm.name
        asrow["auth_alt_id"] = atm.alt_loc or "."
        asrow["auth_comp_id"] = atm.res_name
        asrow["auth_seq_id"] = atm.fragment_id
        asrow["auth_asym_id"] = atm.chain_id
        asrow["type_symbol"] = atm.element
        asrow["Cartn_x"] = atm.position[0]
        asrow["Cartn_y"] = atm.position[1]
        asrow["Cartn_z"] = atm.position[2]
        asrow["occupancy"] = atm.occupancy
        asrow["B_iso_or_equiv"] = atm.temp_factor

        if atm.sig_position:
            asrow["Cartn_x_esd"] = atm.sig_position[0]
            asrow["Cartn_y_esd"] = atm.sig_position[1]
            asrow["Cartn_z_esd"] = atm.sig_position[2]
            asrow["occupancy_esd"] = atm.sig_occupancy
            asrow["B_iso_or_equiv_esd"] = atm.sig_temp_factor

        if atm.U:
            aniso = self.get_table("atom_site_anisotrop")
            anrow = mmCIFRow()
            aniso.append(anrow)

            anrow["id"] = asrow["id"]
            anrow["type_symbol"] = asrow["type_symbol"]
            anrow["label_entity_id"] = asrow["label_entity_id"]
            anrow["pdbx_auth_seq_id"] = asrow["auth_seq_id"]
            anrow["pdbx_auth_comp_id"] = asrow["auth_comp_id"]
            anrow["pdbx_auth_asym_id"] = asrow["auth_asym_id"]
            anrow["pdbx_auth_atom_id"] = asrow["auth_atom_id"]
            anrow["pdbx_auth_alt_id"] = asrow["auth_alt_id"]
            anrow["U[1][1]"] = atm.U[0,0]
            anrow["U[2][2]"] = atm.U[1,1]
            anrow["U[3][3]"] = atm.U[2,2]
            anrow["U[1][2]"] = atm.U[0,1]
            anrow["U[1][3]"] = atm.U[0,2]
            anrow["U[2][3]"] = atm.U[1,2]

            if atm.sig_U:
                anrow["U[1][1]_esd"] = atm.sig_U[0,0]
                anrow["U[2][2]_esd"] = atm.sig_U[1,1]
                anrow["U[3][3]_esd"] = atm.sig_U[2,2]
                anrow["U[1][2]_esd"] = atm.sig_U[0,1]
                anrow["U[1][3]_esd"] = atm.sig_U[0,2]
                anrow["U[2][3]_esd"] = atm.sig_U[1,2]



### <testing>
if __name__ == '__main__':
    import sys

    try:
        path = sys.argv[1]
    except IndexError:
        print "usage: mmCIF.py <mmCIF file path>"
        sys.exit(1)

    cif = mmCIFDictionary()
    cif.load_file(path)
    cif.save_file(sys.stdout)
### </testing>
