package nl.esciencecenter.amuse.distributed.local;

import nl.esciencecenter.amuse.distributed.AmuseConfiguration;
import nl.esciencecenter.amuse.distributed.DistributedAmuseException;
import nl.esciencecenter.octopus.engine.util.StreamForwarder;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import ibis.ipl.server.ServerConnection;

public class Hub {

    private static final long TIMEOUT = 60000;

    private static final Logger logger = LoggerFactory.getLogger(Hub.class);

    private final Process process;

    private final ServerConnection serverConnection;

    private final String address;

    private static String buildCommand(Resource resource, AmuseConfiguration config, String serverAddress) throws DistributedAmuseException {
        StringBuilder command = new StringBuilder();

        //remote host, use ssh
        if (resource.getHostname() != null && !resource.getHostname().equals("localhost")) {
            command.append("ssh ");
            command.append(resource.getUsername());
            command.append("@");
            command.append(resource.getHostname());
            command.append(" -p ");
            command.append(resource.getPort());
            command.append(" ");
        }

        command.append(config.getJava());
        command.append(" -cp ");
        command.append(config.getAmuseHome().getPath());
        command.append("/sandbox/ndrost/distributed/distributed-server.jar");
        command.append(" ibis.ipl.server.Server --remote --hub-only --port 0 --hub-addresses ");
        command.append(serverAddress);

        return command.toString();
    }

    public Hub(Resource resource, AmuseConfiguration config, String serverAddress) throws DistributedAmuseException {
        try {
        String command = buildCommand(resource, config, serverAddress);

        logger.debug("running command " + command);

        process = Runtime.getRuntime().exec(command);
        
        new StreamForwarder(process.getErrorStream(), System.err);

        serverConnection =
                new ServerConnection(process.getInputStream(), process.getOutputStream(), System.out, "Hub at "
                        + resource.getName() + ": ", TIMEOUT, null);

        address = serverConnection.getAddress();
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
