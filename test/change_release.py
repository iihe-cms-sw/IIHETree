import getopt, os, sys

prefix = 'UserCode/IIHETree'
release = 'CMSSW_7_3_0'

argv = sys.argv[1:]
try:
    opts, args = getopt.getopt(argv,"hr:",["release="])
except getopt.GetoptError:
    print 'transform_release.py --release=CMSSW_7_3_0'
    sys.exit(2)
for opt, arg in opts:
    if opt == '--release':
        release = arg

print 'release = %s'%release

directories = ['src','interface']

listOfChangedFiles = []
for dir in directories:
    for filename in os.listdir('%s/%s'%(prefix,dir)):
        changedFile = False
        filepath = '%s/%s/%s'%(prefix,dir,filename)
        if os.path.isfile(filepath):
            file = open(filepath)
            lines_in = file.readlines()
            lines_out = []
            
            for line in lines_in:
                if 'CHOOSE_RELEASE' in line:
                    line = line.replace('\n', '')
                    
                    line = line.replace('/*', '')
                    line = line.replace('*/', '')
                    line = line.replace('//', '')
                    
                    if release in line:
                        changedFile = True
                        line = '// %s'%line
                    else:
                        if 'START' in line:
                            line = '/* %s'%line
                        elif 'END' in line:
                            line = '%s */'%line
                    line = '%s\n'%line
                lines_out.append(line)
                
            file = open(filepath,'w')
            for line in lines_out:
                file.write(line)
        
            if changedFile:
                listOfChangedFiles.append(filepath)
print 'Found that release in %d file(s)'%len(listOfChangedFiles)
for file in listOfChangedFiles:
    print '  %s'%file


