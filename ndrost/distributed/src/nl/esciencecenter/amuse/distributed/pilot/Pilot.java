/*
 * Copyright 2013 Netherlands eScience Center
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package nl.esciencecenter.amuse.distributed.pilot;

import ibis.ipl.Ibis;
import ibis.ipl.IbisFactory;
import ibis.ipl.IbisProperties;

import java.util.Map.Entry;
import java.util.Properties;

import nl.esciencecenter.amuse.distributed.Network;

/**
 * Pilot job. Started when a reservations is made.
 * 
 * @author Niels Drost
 * 
 */
public class Pilot {

    public static void main(String[] arguments) throws Throwable {
        String nodeLabel = "default";
        int slots = 1;
        
        Properties properties = new Properties();
        properties.put(IbisProperties.POOL_NAME, "amuse");
        properties.put(IbisProperties.SERVER_IS_HUB, "false");
        //properties.put("ibis.managementclient", "true");
        //properties.put("ibis.bytescount", "true");

        for(int i = 0; i < arguments.length;i++) {
            if (arguments[i].equalsIgnoreCase("--node-label")) {
                i++;
                nodeLabel = arguments[i];
            } else if (arguments[i].equalsIgnoreCase("--resource-name")) {
                i++;
                properties.put(IbisProperties.LOCATION, arguments[i]);
            } else if (arguments[i].equalsIgnoreCase("--server-address")) {
                i++;
                properties.put(IbisProperties.SERVER_ADDRESS, arguments[i]);
            } else if (arguments[i].equalsIgnoreCase("--hub-addresses")) {
                i++;
                properties.put(IbisProperties.HUB_ADDRESSES, arguments[i]);
            } else if (arguments[i].equalsIgnoreCase("--slots")) {
                i++;
                slots = Integer.parseInt(arguments[i]);
            } else {
                System.err.println("Unknown command line option: " + arguments[i]);
                System.exit(1);
            }
        }
        
        Ibis ibis = IbisFactory.createIbis(Network.IPL_CAPABILITIES, properties, true, null, null, nodeLabel + "," + slots, Network.ONE_TO_ONE_PORT_TYPE);

        System.err.println("running Pilot at location \"" + arguments[0] + "\" for 60 seconds, using properties:");
        for(Entry<Object, Object> entry: properties.entrySet()) {
            System.err.println(entry.getKey() + " = " + entry.getValue());
        }
        
        Thread.sleep(60000);
        
        ibis.end();
    }
}
