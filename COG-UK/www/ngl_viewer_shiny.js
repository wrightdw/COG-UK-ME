document.addEventListener('DOMContentLoaded', function () {
    // create a `stage` object
    var stage = new NGL.Stage('viewportMut', { backgroundColor: 'white' });

    // Create a class with a `atomColor` method that returns a hex color.
    var schemeId = NGL.ColormakerRegistry.addScheme(function (params) {
        this.atomColor = function (atom) {
            return 0xbdb7b7; // default gray
        };
    });

    var tooltip = document.createElement("div");
    Object.assign(tooltip.style, {
        display: "none",
        position: "absolute",
        zIndex: 10,
        pointerEvents: "none",
        backgroundColor: "rgba(0, 0, 0, 0.6)",
        color: "lightgrey",
        padding: "0.5em",
        fontFamily: "sans-serif"
    });
    stage.viewer.container.appendChild(tooltip);

    // Load PDB entry 
    stage.loadFile( "rcsb://6vxx").then(function (o) {
        o.addRepresentation('cartoon', { color: schemeId }); // pass schemeId here
        var pa = o.structure.getPrincipalAxes();
        stage.animationControls.rotate(pa.getRotationQuaternion(), 1500);
    });
    
    //source: https://stackoverflow.com/questions/20798477/how-to-find-index-of-all-occurrences-of-element-in-array
    function getAllIndexes(arr, val) {
        var indexes = [], i;
        for(i = 0; i < arr.length; i++)
            if (arr[i].position === val)
                indexes.push(arr[i].label.substring(1, arr[i].label.length).replace(/[0-9]/g, ''));
        return indexes;
    }

    (function () {
        Shiny.addCustomMessageHandler('vizSpike',
            function (message) {
                var wuhan_data_js = message.a;
                var uploaded_mutations = message.b;
                var uploaded_mutations_positions = message.c;

                var dms_spike_upload_filter_label = [];
                if (uploaded_mutations.length > 0) {
                    dms_spike_upload_filter_label = wuhan_data_js.map(function (item) {
                        return item.position;
                    });
                }

                if (message != '') {
                    var schemeId = NGL.ColormakerRegistry.addScheme(function (params) {
                        this.atomColor = function (atom) {
                            const isPresent = (element) => element == atom.resno;
                            if (Array.isArray(dms_spike_upload_filter_label)) {
                                var indx = dms_spike_upload_filter_label.findIndex(isPresent);
                                if (indx != -1) {
                                  if( message.e[indx].length == 1) {
                                    return message.e.replace('#', '0x');
                                  }
                                  else {
                                    return message.e[indx].replace('#', '0x');
                                  }
                                } else {
                                    return 0xbdb7b7;
                                }
                            } else if (atom.resno == dms_spike_upload_filter_label) {
                                return message.e.replace('#', '0x');
                            } else {
                                return 0xbdb7b7;
                            }
                        };
                    });

                    stage.removeAllComponents();
                    stage.mouseControls.remove("hoverPick");
                    stage.signals.hovered.add(function (pickingProxy) {
                        if (pickingProxy && (pickingProxy.atom || pickingProxy.bond)) {
                            var atom = pickingProxy.atom || pickingProxy.closestBondAtom;
                            var cp = pickingProxy.canvasPosition;
                            var isPresent = (element) => element == atom.resno;
                            var indxFilter = dms_spike_upload_filter_label.findIndex(isPresent);

                            if (indxFilter != -1) {
                                //collect any multiple mutations for a given spike position
                                var indexes = getAllIndexes(wuhan_data_js, atom.resno);

                                tooltip.innerText = "Ref: " + wuhan_data_js[indxFilter].label[0] + "\n" + "Pos: " + wuhan_data_js[indxFilter].label.match(/\d+/) + "\n";
                                if (indexes.length > 1) {
                                    tooltip.innerText += "Alt: " + indexes.join() + "\n";
                                } else {
                                    tooltip.innerText += "Alt: " + indexes[0] + "\n";
                                }

                                switch (message.d) {
                                    case "semantic_score": tooltip.innerText += "Semantic Score: " + wuhan_data_js[indxFilter].semantic_score + "\n";
                                        break;
                                    case "grammaticality": tooltip.innerText += "Grammaticality: " + wuhan_data_js[indxFilter].grammaticality + "\n";
                                        break;
                                    case "evolutionary_index": tooltip.innerText += "Rel. grammaticality: " + wuhan_data_js[indxFilter].evolutionary_index + "\n";
                                        break;
                                    case "entropy": tooltip.innerText += "Entropy: " + wuhan_data_js[indxFilter].entropy + "\n";
                                        break;
                                    case "epitope_max": tooltip.innerText += "Accessibility: " + wuhan_data_js[indxFilter].epitope_max + "\n";
                                        break;
                                }
                                tooltip.innerText += "Evol. Selection: " + wuhan_data_js[indxFilter].evol_selection;
                            } else {
                                if(Array.isArray(uploaded_mutations_positions)){
                                  var indxAllMuts = uploaded_mutations_positions.findIndex(isPresent);
                                  
                                  if (indxAllMuts != -1) {
                                      tooltip.innerText = "Ref: " + uploaded_mutations[indxAllMuts][0] + "\n" + "Pos: " + uploaded_mutations[indxAllMuts].match(/\d+/) + "\n"
                                          + "Alt: " + uploaded_mutations[indxAllMuts].substring(1, uploaded_mutations[indxAllMuts].length).replace(/[0-9]/g, '');
                                  } else {
                                      tooltip.innerText = "Ref: " + "-" + "\n" + "Pos: " + atom.resno + "\n"
                                          + "Alt: " + "-";
                                  }
                                } else {
                                    if (uploaded_mutations_positions == atom.resno) {
                                      tooltip.innerText = "Ref: " + uploaded_mutations[0][0] + "\n" + "Pos: " + uploaded_mutations[indxAllMuts].match(/\d+/) + "\n"
                                          + "Alt: " + uploaded_mutations[0].substring(1, uploaded_mutations[0].length).replace(/[0-9]/g, '');
                                    } else {
                                        tooltip.innerText = "Ref: " + "-" + "\n" + "Pos: " + atom.resno + "\n"
                                            + "Alt: " + "-";
                                    }
                                }
                            }

                            tooltip.style.bottom = cp.y + 3 + "px";
                            tooltip.style.left = cp.x + 3 + "px";
                            tooltip.style.display = "block";
                        } else {
                            tooltip.style.display = "none";
                        }
                    });


                    stage.loadFile("rcsb://6vxx").then(function (o) {
                        o.addRepresentation(message.f, { color: schemeId }); // pass schemeId here
                        o.autoView();
                    });
                }
            });


    })();

}
);