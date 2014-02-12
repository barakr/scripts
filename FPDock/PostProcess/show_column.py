import re
import sys

def usage_and_exit():
    print "\t", sys.argv[0], "<fname> <column_name1> [column_name2] ..."
    exit()

def get_column_id(header, names):
    """ returns the column ids of each name in names, within header """
#    print names
    cols = header.split()
    retval=[]
    for name in names:
#        print name,
        is_found = False
        for id,col in enumerate(cols):
            if(col == name):
                retval.append(id)
                is_found = True
        if(not is_found):
            raise KeyError("column " + name + " not found")
#    print
    return retval

def get_column_id_dict(header, names):
    """ returns the column ids of each name in names, within header
    -1 for all columns that are not found
    """
#    print names
    cols = header.split()
    dict={}
    for name in names:
#        print name,
        is_found = False
        for id,col in enumerate(cols):
            if(col == name):
                dict[name] = id
                is_found = True
        if(not is_found):
            dict[name] = -1
#    print
    return dict

def read_scorefile(fname, captions, entries):
    '''
    Read the columns 'captions' from whitespace delimeted file 'fname'.
    The first line in the file is assumed to be a header with caption names.

    fname - the file name
    captions - the captions of column to be read
    entries[in/out] - each entry mapping from caption
    to value, is appended into the entries list.
    '''
    F=open(fname,'r')
    header=F.readline()
    caption2id = get_column_id_dict(header, captions)
#    print fname, caption2id
    for line in F:
        line_as_list = line.split()
        entry={'fname' : fname}
        try:
            for caption in captions:
                col_id = caption2id[caption]
                if col_id > -1:
                    value = line_as_list[col_id]
                else:
                    value = "NA"
                entry[caption] = value
            entries.append(entry)
        except:
            print >>sys.stderr, "WARNING: Skipping line: '", line.rstrip('\n') , "'"


def main():
    if(len(sys.argv) <= 2):
        usage_and_exit()
    F=open(sys.argv[1],'r')
    header=F.readline()
    ids = get_column_id(header, sys.argv[2:])
    for line in F:
        data = line.split()
        for id in ids:
            try:
                print (data[id]),
            except:
                print "ERR",
        print

if (__name__ == "__main__"):
    main()
