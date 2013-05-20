package nl.esciencecenter.amuse.distributed.local;

import java.net.URI;
import java.util.List;

import nl.esciencecenter.amuse.distributed.AmuseConfiguration;
import nl.esciencecenter.amuse.distributed.DistributedAmuseException;
import nl.esciencecenter.octopus.Octopus;
import nl.esciencecenter.octopus.credentials.Credential;
import nl.esciencecenter.octopus.engine.util.StreamForwarder;
import nl.esciencecenter.octopus.jobs.Job;
import nl.esciencecenter.octopus.jobs.JobDescription;
import nl.esciencecenter.octopus.jobs.Scheduler;
import nl.esciencecenter.octopus.util.JavaJobDescription;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import ibis.ipl.server.ServerConnection;

public class Hub {

    private static final long TIMEOUT = 60000;

    private static final Logger logger = LoggerFactory.getLogger(Hub.class);

    private final Job job;

    private final ServerConnection serverConnection;

    private final String address;

    private static JavaJobDescription createJobDesciption(Resource resource, String[] hubAddresses)
            throws DistributedAmuseException {
        JavaJobDescription result = new JavaJobDescription();

        result.setInteractive(true);

        AmuseConfiguration configuration = resource.getConfiguration();
        
        result.setExecutable(configuration.getJava());

        //classpath
        List<String> classpath = result.getJavaClasspath();
        classpath.add(configuration.getAmuseHome().getPath() + "/sandbox/ndrost/distributed/distributed-server.jar");
        classpath.add(configuration.getAmuseHome().getPath() + "/sandbox/ndrost/distributed");

        result.setJavaMain("ibis.ipl.server.Server");

        //arguments
        List<String> javaArguments = result.getJavaArguments();
        javaArguments.add("--remote");
        javaArguments.add("--hub-only");
        javaArguments.add("--port");
        javaArguments.add("0");

        String hubs = null;
        for (String hub : hubAddresses) {
            if (hubs == null) {
                hubs = hub;
            } else {
                hubs = hubs + "," + hub;
            }
        }

        if (hubs != null) {
            javaArguments.add("--hub-addresses");
            javaArguments.add(hubs);
        }

        return result;
    }

    public Hub(Resource resource, AmuseConfiguration config, String[] hubs, Octopus octopus) throws DistributedAmuseException {
        try {
            JobDescription jobDescription = createJobDesciption(resource, hubs);

            logger.debug("starting hub with job description " + jobDescription + " with arguments " + jobDescription.getArguments());

            Scheduler scheduler;

            if (resource.isLocal()) {
                scheduler = octopus.jobs().getLocalScheduler();
            } else {
                //FIXME, replace with octopus.credentials.getDefaultCredential("ssh") once this is implemented
                String username = System.getProperty("user.name");
                Credential credential =
                        octopus.credentials().newCertificateCredential("ssh", null, "/home/" + username + "/.ssh/id_rsa",
                                "/home/" + username + "/.ssh/id_rsa.pub", username, "");

                URI uri = new URI("ssh", resource.getUsername(), resource.getHostname(), resource.getPort(), null, null, null);
                scheduler = octopus.jobs().newScheduler(uri, credential, null);
                
                logger.debug("starting hub using scheduler " + scheduler);
            }

            job = octopus.jobs().submitJob(scheduler, jobDescription);
            
            logger.debug("started job "  + job);
            
            while(job.getStdout() == null) {
                //FIXME: workaround for interactive jobs that do not start immediately
                logger.warn("waiting until interactive job has started (workaround)");
                Thread.sleep(100);
            }

            new StreamForwarder(job.getStderr(), System.err);

            serverConnection =
                    new ServerConnection(job.getStdout(), job.getStdin(), System.out, "Hub at " + resource.getName() + ": ",
                            TIMEOUT, null);

            address = serverConnection.getAddress();
            
            logger.debug("hub on " + resource.getName() + " has address " + address);
        } catch (Exception e) {
            throw new DistributedAmuseException("cannot start hub on " + resource.getName(), e);
        }
    }

    String getAddress() {
        return address;
    }

    void stop() {
        serverConnection.closeConnection();
    }

}
