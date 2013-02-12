package nl.esciencecenter.estars;

import java.io.IOException;
import java.lang.System;
import java.net.InetSocketAddress;
import java.net.Socket;
import java.net.SocketAddress;
import java.nio.channels.SocketChannel;

//private class AmuseMessage{} .....

//generated Worker class
class Worker {

    private final AmuseMessage request;
    private final AmuseMessage reply;
    
    
    
    Worker() {
        this.request = new AmuseMessage();
        this.reply = new AmuseMessage();
    }

    
    //generated
    private void handleCall() {
        if(request.getFunctionID() == 1060350414) {
            reply.setHeader(1, 0, 0, 0, 0, 0);
            System.err.println("echoing int");
            
        }
        
    }

    
    private void runSockets() {
        
    }
    
    
    
    
    public static void main(String[] arguments) throws IOException {
        SocketChannel channel;

        System.err.println("Java eStars Worker");
        for (String argument : arguments) {
            System.err.println("argument: " + argument);
        }
        
        if (arguments.length == 0) {
            System.err.println("No arguments to worker. expected a socket port number");
            System.exit(1);
        }
        
        int port = Integer.parseInt(arguments[0]);
        
        channel = SocketChannel.open(new InetSocketAddress(port));
        
        while(true) {
            AmuseMessage request = new AmuseMessage();
            request.readFrom(channel);
            
            System.err.println("got message " + request.toContentString());
            
            AmuseMessage reply = new AmuseMessage();
            
            if(request.getFunctionID() == 1060350414) {
                reply.setCallID(request.getCallID());
                reply.setFunctionID(request.getFunctionID());
                reply.setCallCount(request.getCount());
                
            }
            
            
            request.clear();
            
            
            
            
            request.writeTo(channel);
        }

    }
}
