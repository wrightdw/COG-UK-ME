function _get_spike(data) {
    return (
        _.chain(data).toPairs().filter((d) => d[1].G === 'S').fromPairs().value()
    )
}

function _filter_baa(data) {
    return _.chain(data).toPairs().filter((d) => 'baa' in d[1]).fromPairs().value();
}

function _process_summary_json(data) {
    var sc2_seq = data['SC2'];
    delete data['SC2'];
    data = _.chain(data).toPairs().filter(d => d[1].busted).fromPairs().value();
    data = [sc2_seq, data][1];

    return data;
}

function _species(summary_json) {
    return (
        _.sortBy(_.filter(_.uniq(_.map(summary_json[_.first(_.keys(summary_json))]['tree_tags'], (d) => d[0])), (d) => d.length))
    )
}

function _prime_annotation(summary_json) {
    return (
        _.map(summary_json[_.first(_.keys(summary_json))]["prime-properties"], (d) => d.split(" ")[2])
    )
}

function _site2segment(summary_json) {
    let r = {};
    _.each(summary_json, (d, s) => {
        _.each(d.map, (c) => {
            if (c >= 0) {
                r[c] = s;
            }
        });
    });

    return r;
}

function _site2alignment(annotation_json, site2segment, summary_json) {
    let r = {};
    let offset = 1;
    let last_segment = "";
    let only_sites = _.sortBy(_.map(annotation_json, (d, k) => [k, d.index]), (s) => +s.k);
    _.each(only_sites, (d) => {
        if (site2segment[d[0]] != last_segment) {
            if (last_segment && last_segment.length) {
                offset += summary_json[last_segment].map.length;
            }
            last_segment = site2segment[d[0]];
        }
        r[d[0]] = offset + d[1];
    });

    return r;
}

function _conserved_sites(annotation_json, species) {
    return (
        _.filter(_.map(annotation_json, (v, k) => [k, v]), (s) => _.uniq(_.flatten(_.map(species, (sp) => _.uniq(_.keys(s[1].baa[sp]))))).length == 1, (d) => [+d[0]])
    )
}

function _selected_sites(annotation_json, site2alignment, prime_annotation) {
    return (
        _.map(_.filter(_.map(annotation_json, (v, k) => [k, v]), (s) => s[1].bMEME.p <= 0.05), (d) => { return { 'coordinate': +d[0], 'p': d[1].bMEME.p, 'meme_weight': d[1].bMEME["w+"], 'gene': d[1].G + "/" + d[1].S, 'site': site2alignment[d[0]], 'fel_p': d[1].bFEL.p, 'fel_n': d[1].bFEL.b < d[1].bFEL.a, 'prime': _.map(_.filter(_.map(d[1].bPRIME ? d[1].bPRIME.p : [], (p, i) => [i, p, i > 0 ? d[1].bPRIME['lambda'][i - 1] : null]), (pr) => pr[1] < 0.05), (d) => { if (d[0] > 0) d[0] = prime_annotation[d[0] - 1]; else d[0] = "Overall"; return d; }) }; })
    )
}

function _negative_sites(annotation_json, species, site2alignment) {
    return (
        _.map(_.filter(_.map(annotation_json, (v, k) => [k, v]), (s) => _.size(s[1].bcdn[species[1]]) > 1 && s[1].bFEL.p <= 0.05 && s[1].bFEL.a > s[1].bFEL.b), (d) => [+d[0], d[1].bFEL.b / d[1].bFEL.a, d[1].bFEL.p, site2alignment[d[0]]])
    )
}

function _credibility_heatmap(annotation_json) {
    return _.map(_.filter(_.map(annotation_json, (d, s) => {
        if ("evo" in d) {
            return [s,
                d3.sum(_.map(d["evo"], (d) => -d * Math.log2(d)))

            ];
        }
        return [s, null];
    }), (d) => d[1]), (d) => { return { 'coordinate': d[0], 'entropy': d[1] } });
}

function _convert_to_spike_coordinates(res_obj, mode) {
    //var indxs = [];
    if (mode == 'default') {
        for (var i = 0; i < res_obj.length; i++) {
            res_obj[i][0] = (+[res_obj[i][0]] - 21562) + 1;
        }
    } else if (mode == 'coordinate') {
        for (var i = 0; i < res_obj.length; i++) {
            res_obj[i].coordinate = (+[res_obj[i].coordinate] - 21562) + 1;
        }
    }
    return res_obj;
}

function _get_indexes(res_obj, mode) {
    var indxs = [];
    if (res_obj != undefined) {
        if (mode == 'default') {
            for (var i = 0; i < res_obj.length; i++) {
                indxs.push(+[res_obj[i][0]]);
            }
        } else if (mode == 'coordinate') {
            for (var i = 0; i < res_obj.length; i++) {
                indxs.push(+[res_obj[i].coordinate]);
            }
        }
    }

    return indxs;
}

function _get_amino_acid_locations(nucl_indxs) {
    var amino_indxs = [];
    for (var i = 0; i < nucl_indxs.length; i++) {
        if (nucl_indxs[i] <= 3) amino_indxs.push(1);
        else {
            var num = ~~(nucl_indxs[i] / 3);
            if (num === 0) {
                amino_indxs.push(num);
            } else {
                amino_indxs.push(num + 1);
            }
        }
    }

    return amino_indxs;
}

function _assign_selection_to_sites(objArray, sites_to_assign, evol_sel) {
    var array = typeof objArray != 'object' && objArray !== null ? JSON.parse(objArray) : objArray;
    var str = '';
    var line = '';
    for (var col in array[0]) line += col + ','
    line += evol_sel; //add new header for evolutionary selection column
    str += line + '\n';

    if (sites_to_assign.length > 0) {
        sites_to_assign.sort(function (a, b) { return a - b });
        sites_to_assign = _get_amino_acid_locations(sites_to_assign);

        for (var i = 0; i < array.length; i++) {
            line = '';
            for (var j in array[i]) {
                line += array[i][j] + ',';
            }
            str += line;

            line = '';
            const isPresent = (element) => element == array[i].position; //locate position column
            var indxPos = sites_to_assign.findIndex(isPresent);
            if (indxPos != -1) {
                line += 'Yes'; //assign value to the last column added above
            } else {
                line += 'No';
            }
            if (i < array.length - 1) str += line + '\n';
            else str += line;
        }
    } else {
        for (var i = 1; i < array.length; i++) {
            str += array[i];
            line = ',' + 'No';
            if (i < array.length - 1) str += line + '\n';
            else str += line;
        }
    }

    return str;
}

function _assign_metadata_to_sites_csv(csvArray, sites_to_assign, credibility_heatmap, tag, metadata) {
    var str = '';
    var line = '';

    for (var col in csvArray[0]) line += csvArray[0][col] + ','
    line += tag; //add new header for evolutionary selection column
    str += line + '\n';

    if (sites_to_assign.length > 0) {
        sites_to_assign.sort(function (a, b) { return a - b });
        sites_to_assign = _get_amino_acid_locations(sites_to_assign);

        for (var i = 1; i < csvArray.length; i++) {
            line = '';
            for (var j in csvArray[i]) {
                line += csvArray[i][j] + ',';
            }
            str += line;

            line = '';
            const isPresent = (element) => element == csvArray[i][3]; //locate position column
            var indxPos = sites_to_assign.findIndex(isPresent);
            if(indxPos != -1) {
                switch (metadata) {
                    case 'selection':
                        line += 'Yes'; //assign value to the last column added above
                        break;
                    case 'entropy':
                        line += credibility_heatmap[indxPos].entropy;
                        break;
                }
            } else {
                switch (metadata) {
                    case 'selection':
                        line += 'No'; //assign value to the last column added above
                        break;
                    case 'entropy':
                        line += "NA";
                        break;
                }
            }
            if(i < csvArray.length - 1) str += line + '\n';
            else str += line;
        }
    } else {
        for (var i = 1; i < csvArray.length; i++) {
            str += csvArray[i];
            switch (metadata) {
                    case 'selection':
                        line = ',' + 'No'; //assign value to the last column added above
                        break;
                    case 'entropy':
                        line += ',' + "NA";
                        break;
                }
            if(i < csvArray.length - 1) str += line + '\n';
            else str += line;
        }
    }

    return str;
}

function _main(annotation_json, summary_json, wuhan_data) {
    annotation_json = _filter_baa(annotation_json);
    annotation_json = _get_spike(annotation_json);
    summary_json = _process_summary_json(summary_json);

    var species = _species(summary_json);
    var prime_annotation = _prime_annotation(summary_json);
    var site2segment = _site2segment(summary_json);
    var site2alignment = _site2alignment(annotation_json, site2segment, summary_json);

    var conserved_sites = _conserved_sites(annotation_json, species);
    var selected_sites = _selected_sites(annotation_json, site2alignment, prime_annotation);
    var negative_sites = _negative_sites(annotation_json, species, site2alignment);
    var credibility_heatmap = _credibility_heatmap(annotation_json);

    credibility_heatmap = _convert_to_spike_coordinates(credibility_heatmap, mode = 'coordinate');
    conserved_sites = _convert_to_spike_coordinates(conserved_sites, mode = 'default');
    negative_sites = _convert_to_spike_coordinates(negative_sites, mode = 'default');
    positive_sites = _convert_to_spike_coordinates(selected_sites, mode = 'coordinate');

    var dms_spike_upload_filter = wuhan_data; // used when we want to annotate the full data (all mutations per aa position) and extract them into a CSV file
    var dms_spike_upload_filter_csv = [];
    if (dms_spike_upload_filter.length > 0) {
        dms_spike_upload_filter_csv = _assign_selection_to_sites(dms_spike_upload_filter, _get_indexes(conserved_sites, 'default'), "conserved");

        const csv_to_array = (data, delimiter = ',', omitFirstRow = false) =>
            data
                .slice(omitFirstRow ? data.indexOf('\n') + 1 : 0)
                .split('\n')
                .map(v => v.split(delimiter));

        dms_spike_upload_filter_csv = _assign_metadata_to_sites_csv(csv_to_array(dms_spike_upload_filter_csv), _get_indexes(positive_sites, 'coordinate'), null, tag = "positive", metadata = "selection");

        dms_spike_upload_filter_csv = _assign_metadata_to_sites_csv(csv_to_array(dms_spike_upload_filter_csv), _get_indexes(negative_sites, 'default'), null, tag = "negative", metadata = "selection");

        dms_spike_upload_filter_csv = _assign_metadata_to_sites_csv(csv_to_array(dms_spike_upload_filter_csv), _get_indexes(credibility_heatmap, 'coordinate'), credibility_heatmap, tag = "entropy", metadata = "entropy");
    }

    return dms_spike_upload_filter_csv;
}