#!/usr/bin/python
# SEE SECURITY HAZARD WARNING BELOW
import show_column
import sys
import subprocess
import re
import getopt

def get_cmdline_options():
    def usage_and_exit():
        print sys.argv[0], '[-n <number_of_tops> -e <True/False> -c <captions> -s <sort_col>] file1 [file2 ...]'
        print "captions - colons delimited extra column captions to display, e.g. 'rmsBB:rmsALL'"
        print "sort_col - column by which to sort"
        sys.exit(2)
    if len(sys.argv)==1:
        usage_and_exit()
#    fnames=sys.argv[1:]
    ntop=10
    is_extract=False
    captions=set(['reweighted_sc', 'description', 'rmsBB', 'bestRMS_4mer_all'])
    sort_col="reweighted_sc"
    try:
        opts,fnames = getopt.getopt(sys.argv[1:],"h:n:e:c:s:",["ntop=","extract=","captions=","sort="])
        for opt, arg in opts:
            if opt == ('-h', "--help"):
                usage_and_exit()
            elif opt in ("-n", "--ntop"):
                ntop = int(arg)
            elif opt in ("-e", "--extract"):
                is_extract = arg.lower() in ["1", "true", "t", "yes", "y"]
            elif opt in ("-c", "--captions"):
                captions = captions | set(arg.split(":"))
            elif opt in ("-s", "--sort_col"):
                sort_col = arg
                captions = captions | set([sort_col]) # make sure it is there
    except getopt.GetoptError:
        usage_and_exit()
    return ntop, is_extract, captions, sort_col, fnames

ntop, is_extract, captions,sort_col,fnames = get_cmdline_options()
entries = []
for fname in fnames:
    try: show_column.read_scorefile(fname, captions,entries)
    except:
        print "Error reading %s - skipping" % fname
entries.sort(key = lambda entry: float(entry["reweighted_sc"]))
fraction = len(entries)/3
entries = entries[0:fraction] # filter out high energy ones
for i,e in enumerate(entries):
    e["rank"]=i
#    e.append(i)
entries.sort(key = lambda entry: float(entry[sort_col]))
for caption in entries[0].keys(): print caption,
print
for entry in entries[0:ntop]:
    for val in entry.values(): print val,
    print
    if(is_extract):
        score_fname = entry['fname']
        silent_fname = re.sub(r'(.*)score_(.*)\.sc$',r'\1decoys_\2.silent', score_fname)
        # NOTE: next line is a SECURITY HAZARD due to possibly malicious content of entry[fname] or entry[description]
        #       do not use unless by a trusted user
        subprocess.call("/netapp/sali/barak/Rosetta/rosetta_source/bin/extract_pdbs.default.linuxgccrelease " +
                        "-database /netapp/sali/barak/Rosetta/rosetta_database/ -silent_read_through_errors -in:file:silent "
                        + silent_fname + " -tags " + entry['description'], shell=True)
        subprocess.call("/bin/mv %s.pdb %s.%s.pdb" %
                        (entry['description'],
                         entry['description'],
                         silent_fname),
                        shell=True)
