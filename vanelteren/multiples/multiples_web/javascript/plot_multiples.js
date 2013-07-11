/**
 *
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the 'License'); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an 'AS IS' BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
$(document).ready(function(){

  var client, destination;

    var url = 'ws://localhost:61623/stomp';
    var login = 'guest';
    var passcode = 'password';
    destination = '/queue/data';
    
    var plot = $.plot(
        "#messages", 
        [ 
            {
                data: [[0,0]],
                points: { show: true }
            } 
        ], 
        {
                series: {
                    shadowSize: 0	// Drawing is faster without shadows
                },
                yaxis: {
                    min: -2,
                    max: 2
                },
                xaxis: {
                    min: -2,
                    max: 2
                }
        }
    );
    
    function plot_encounter(data1, data2, jqElement) {
        var points = [];
        var points2 = [];
        var maxX = -1e20;
        var minX = 1e20;
        var maxY = -1e20;
        var minY = 1e20;
        for (var i = 0; i < data1.length; ++i) {
            var element = data1[i];
            
            var x = element.x;
            var y = element.y;
            var vx = 0.001 * element.vx;
            var vy = 0.001 * element.vy;
            
            meanX += x;
            meanY += y;
            
            points.push([x,y,vx,vy]);
            
            minX = Math.min(minX, x);
            maxX = Math.max(maxX, x);
            
            minY = Math.min(minY, y);
            maxY = Math.max(maxY, y);
            
            minX = Math.min(minX, x + vx);
            maxX = Math.max(maxX, x + vx);
            
            minY = Math.min(minY, y + vy);
            maxY = Math.max(maxY, y + vy);
            
        }
        
        for (var i = 0; i < data2.length; ++i) {
            var element = data2[i];
            
            var x = element.x;
            var y = element.y;
            var vx = 0.001 * element.vx;
            var vy = 0.001 * element.vy;
            
            meanX += x;
            meanY += y;
            
            points2.push([x,y,vx,vy]);
            
            minX = Math.min(minX, x);
            maxX = Math.max(maxX, x);
            
            minY = Math.min(minY, y);
            maxY = Math.max(maxY, y);
            
            minX = Math.min(minX, x + vx);
            maxX = Math.max(maxX, x + vx);
            
            minY = Math.min(minY, y + vy);
            maxY = Math.max(maxY, y + vy);
            
        }
        
        var meanX = (maxX - minX) / 2.0;
        var meanY = (maxY - minY) / 2.0;
        
        $.plot(
            jqElement, 
            [ 
                {
                    data: points,
                    points: { show: true }
                } ,
                {
                    data: points2,
                    points: { show: true }
                } 
            ], 
            {
                series: {
                    shadowSize: 0	// Drawing is faster without shadows
                    ,
                    direction: {
                        show: true
                    }
                },
                xaxis: {
                    min: minX,
                    max: maxX
                },
                yaxis: {
                    min: minY,
                    max: maxY
                }
            }   
        );
    }

    var queue = [];
    setInterval(function() {
        if(queue.length == 0) {
            return;
        }
        var message = queue.shift();
        if(message.type == 'particles') {
            
            $('#title').text(message['time-str']);
            $('#nmultiples').text(message['n-multiples']);
            
            var points = [];
            for (var i = 0; i < message.x.length; ++i) {
                points.push([message.x[i],message.y[i]]);
            }
            var points2 = [];
            for (var i = 0; i < message.multiples.x.length; ++i) {
                points2.push([message.multiples.x[i],message.multiples.y[i]]);
            }
            plot.setData([{data:points, points:{show:true}}, {data:points2, points:{show:true}}]);
            plot.draw()
        } else if (message.type == 'encounter'){
            var table = $('#encounters');
            var row = $('<tr></tr>');
            var before = $('<td></td>');
            before.addClass('encounter-view');
            row.append(before);
            //var after = $('<td></td>');
            //after.addClass('encounter-view');
            //row.append(after);
            var time = $('<td></td>');
            time.text(message['time-str']);
            row.append(time);
            table.prepend(row);
            plot_encounter(message.before, message.after, before);
            
            if(table.children('tbody').children().size() > 20) {
                table.children('tbody').children().last().remove();
            }
        }
    }, 50);
    
    client = Stomp.client(url); //, ['binary']);

    // this allows to display debug logs directly on the web page
    client.debug = function(str) {
      //$("#debug").append(document.createTextNode(str + "\n"));
      //console.log(str);
    };
    // the client is notified when it is connected to the server.
    var onconnect = function(frame) {
        client.debug("connected to Stomp");
        client.subscribe(
            destination, 
            function(message) {
                var message = JSON.parse(message.body);
                queue.push(message);
            /*if(message['time-value'] > 6.5 && message['time-value'] < 7.35) {
                
            }*/
            }
            ,
            {
                'browser':true,
                'browser-end':false,
                'include-seq':'seq',
                'from-seq':0
            }
        );
    };
    client.connect(login, passcode, onconnect);


});
