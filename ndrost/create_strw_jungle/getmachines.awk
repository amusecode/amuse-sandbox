# getscreen all | awk -f getmachines.awk
BEGIN{}

!(/Windows/ || /Virtual/ || /virthost/ || /backblaze/ || /gpu-3/ || /opensim1/ || /para/) {
 platform="";line="";
 ssh="ssh -o StrictHostKeyChecking=no -o BatchMode=yes "$1" sh"
 print " uname -i" |&ssh ;
 ssh |& getline platform;
 if (platform == "x86_64" ) {
   print " java -version |& head -n 2" |& ssh;
   ssh |& getline line;
   ssh |& getline line;
   isopenjdk=match(line,/OpenJDK/) || match(line,/SE Runtime/);
   if (isopenjdk==1) {
     print "source $MODULESHOME/init/sh" |& ssh
     print "module load mpich2-x86_64" |& ssh
     print "if which mpirun > /dev/null; then echo exists ;\
     else echo does not exist; fi" |& ssh;
     ssh |& getline line;
     gotmpirun=0;if(line=="exists") gotmpirun=1;
     if(gotmpirun==1) {
       print "cat /proc/meminfo | grep MemTotal | awk '{print $2}'" |& ssh;
       ssh |& getline line; 
       mem=line;
       print "cat /proc/cpuinfo | grep processor | wc -l" |& ssh;
       ssh |& getline line; 
       cpu=line;
      # print platform, isopenjdk,gotmpirun
       print $1,mem,cpu;
       }
     }
   }
 close(ssh);
 }
 
 
END{}
