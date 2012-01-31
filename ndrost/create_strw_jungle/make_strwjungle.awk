BEGIN{
print "\
default.file.adaptors = SftpTrilead,SshTrilead,CommandlineSSH,Local\n\
default.java.path = java\n\
#default.user.name = pelupes\n\
default.latitude = 52.1686\n\
default.longitude = 4.4598\n\
default.amuse.home= /disks/koppoel2/ndrost/amuse-svn\n\
default.mpirun= /disks/koppoel2/ndrost/amuse-svn/lib/ibis/strw.mpd.wrapper.script\n\
"

}

!/^#/ {hostname=$1; 
print hostname".support.uri = ssh://"hostname".strw.leidenuniv.nl";
print hostname".job.uri = ssh://"hostname".strw.leidenuniv.nl";
print hostname".cores = "$3
print hostname".memory = "int($2/1024)
print hostname".latitude = "52.16873+0.005*(rand()-0.5)
print hostname".longitude = "4.45811+0.005*(rand()-0.5)
}

END{}
