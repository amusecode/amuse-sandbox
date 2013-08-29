package nl.esciencecenter.amuse.distributed.resources;

import ibis.ipl.server.ServerConnection;

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
import nl.esciencecenter.octopus.jobs.Streams;
import nl.esciencecenter.octopus.util.JavaJobDescription;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

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
        //javaArguments.add("--hub-only");
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

            logger.debug("starting hub with job description " + jobDescription + " with arguments "
                    + jobDescription.getArguments());

            Scheduler scheduler;

            if (resource.isLocal()) {
                scheduler = octopus.jobs().getLocalScheduler();
            } else {
                Credential credential = octopus.credentials().getDefaultCredential("ssh");

                URI uri = new URI("ssh", resource.getUsername(), resource.getHostname(), resource.getPort(), null, null, null);
                scheduler = octopus.jobs().newScheduler(uri, credential, null);

                logger.debug("starting hub using scheduler " + scheduler);
            }

            job = octopus.jobs().submitJob(scheduler, jobDescription);

            logger.debug("started job " + job);

            Streams streams = octopus.jobs().getStreams(job);

            new StreamForwarder(streams.getStderr(), System.err);

            serverConnection =
                    new ServerConnection(streams.getStdout(), streams.getStdin(), System.out, "Hub at " + resource.getName()
                            + ": ", TIMEOUT, null);

            address = serverConnection.getAddress();

            logger.debug("hub on " + resource.getName() + " has address " + address);
        } catch (Exception e) {
            throw new DistributedAmuseException("cannot start hub on " + resource.getName(), e);
        }
    }

    public String getAddress() {
        return address;
    }

    void stop() {
        serverConnection.closeConnection();
    }

}
