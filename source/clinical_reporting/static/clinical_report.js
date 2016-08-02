$(function() {
    //$("[rel=tooltip]").tooltip({animation: false});

    $(".tooltip-link").tooltip({ show: { effect: "blind", duration: 0 } });

    if (msieversion() == 0) {
        //$('table.tableSorter').tableSort();

    } else {
        //$('.coverage_plot').hide();
    }

    //if ($.browser.msie) {
    //
    //}
});

var table_short_clones = [];
var table_full_clones = [];
var key_or_target = null;

var minAF = 0;

var showBlacklisted = false;

var mutPlotInfo;

var parameters = [];

$(function() {
    var tables_short = $('.table_short');
    var tables_full = $('.table_full');
    for (var t = 0; t < tables_short.length; t++){
        var table_short = $(tables_short[t]);
        var table_full = $(tables_full[t]);
        full_table_parent = table_full.parent();
        if (full_table_parent.length > 0) {
            var switch_id = full_table_parent[0].id.split('_')[0];
            var switch_el = $('#' + switch_id + '_switch');
            if (switch_id == 'seq2c' && switch_el[0]) key_or_target = switch_el.html().indexOf('target') != -1 ? 'target' : 'key';
            else if (switch_id == 'variants') write_to_excel(table_full);

            table_short_clones.push({'id_': switch_id, 'table': $(table_short).clone()});
            table_full_clones.push({'id_': switch_id, 'table': $(table_full).clone()});

            table_short.remove();
            table_full.remove();

            if (table_short.find('tr').length > 0) {
                reduceClick('extend_link_' + switch_id);
            } else {
                extendClick('extend_link_' + switch_id);
            }
        }
    }
    if ($('#mut_af_textbox')[0]) filterMutationsByAF($('#mut_af_textbox')[0].value);

    $('#variants_table_controls').width($('#report_table_mutations').width() - 5);
    $('#download_mut_table').show();
    createSelectVariants();
    //if (msieversion() == 0) {
    //    $('table.tableSorter.table_short').tableSort();
    //}
});

function extendClick(switch_id) {
    if (switch_id[0].id) switch_id = switch_id[0].id;
    showBlacklisted = switch_id.search('incidentalome') != -1;
    // Showing full
    switch_id = switch_id.split("_");
    var table_id = switch_id[switch_id.length - 1];
    var switch_el = $('#' + table_id + '_switch');
    if (table_id == 'seq2c') {
      switch_el.html('<a class="dotted-link" id="reduce_link_' + table_id + '" onclick="reduceClick($(this))">' + key_or_target + ' genes</a> / <span>all genes</span>')
    } else {
        switchElContent = '<a class="dotted-link" id="reduce_link_' + table_id + '" onclick="reduceClick($(this))">known, likely</a> / ';
        if (showBlacklisted) {
            switchElContent += '<a class="dotted-link" id="extend_link_' + table_id + '" onclick="extendClick($(this))">+ unknown</a>';
            //switchElContent += '<span id="incidentalome_span">incidentalome</span>';
        }
        else {
            switchElContent += '<span>+ unknown</span>';
            //switchElContent += '<a class="dotted-link" id="extend_link_incidentalome_' + table_id + '" onclick="extendClick($(this))">incidentalome</a>';
        }
        switch_el.html(switchElContent)
    }
    var table_div = $('#' + table_id + '_table_div');
    var table_short = table_div.find('.table_short');
    if (table_short) table_short.remove();
    var table_full = table_div.find('.table_full');
    if (!table_full[0]) {
        var table_full_clone = null;
        for (var t = 0; t < table_full_clones.length; t++) {
          if (table_full_clones[t].id_ == table_id) {
            table_full_clone = table_full_clones[t].table;
            table_full_clones.splice(t, 1);
          }
        }
        table_div.prepend(table_full_clone);
        table_full = table_div.find('.table_full');
    }
    if (table_full) {
        if (table_id == 'variants') {
            $(table_full).css('height', '');
            $('.table_full#report_table_mutations tr').each(function() {
                checkBlacklisted(this, showBlacklisted);
                checkAF(this, minAF);
                if ($('#show_variants_all')[0]) {
                    showVariantsBySensitivity('show_variants_all');
                    showVariantsByType('show_variants_all');
                }
            });
        }
        $(table_full).show();
        if (msieversion() == 0) {
            table_full.tableSort();
        }
    }
    if (table_full_clone) table_full_clones.push({'id_': table_id, 'table': table_full_clone.clone()});
    if (mutPlotInfo) showPlotWithInfo(mutPlotInfo, minAF);
}

function reduceClick(switch_id) {
    if (switch_id[0].id) switch_id = switch_id[0].id;
    showBlacklisted = false;
    switch_id = switch_id.split("_");
    var table_id = switch_id[switch_id.length - 1];
    // Showing reduced
    var switch_el = $('#' + table_id + '_switch');
    if (table_id == 'seq2c') {
      switch_el.html('<span>'  + key_or_target + ' genes</span> / <a class="dotted-link" id="extend_link_' + table_id + '" onclick="extendClick($(this))">all genes</a>')
    }
    else {
        html = '<span>known, likely</span> / <a class="dotted-link" id="extend_link_' + table_id + '" ' +
          'onclick="extendClick($(this))">+ unknown</a>'
        //html += '<a class="dotted-link" id="extend_link_incidentalome_' + table_id + '" onclick="extendClick($(this))">incidentalome</a>'
      switch_el.html(html)
    }
    var table_div = $('#' + table_id + '_table_div');
    var table_full = table_div.find('.table_full');
    if (table_full) table_full.remove();
    var table_short_clone = null;
    for (var t = 0; t < table_short_clones.length; t++) {
      if (table_short_clones[t].id_ == table_id) {
        table_short_clone = table_short_clones[t].table;
        table_short_clones.splice(t, 1);
      }
    }
    table_div.prepend(table_short_clone);
    var table_short = table_div.find('.table_short');
    if (table_short) {
      if (table_id == 'variants') {
        $(table_short).css('height', '');
        $('.table_short#report_table_mutations tr').each(function() {
            checkAF(this, minAF);
            if ($('#show_variants_all')[0]) {
                showVariantsBySensitivity('show_variants_all');
                showVariantsByType('show_variants_all');
            }
        });
      }
      table_short.show();
      if (msieversion() == 0) {
          table_short.tableSort();
      }
    }
    table_short_clones.push({'id_': table_id, 'table': table_short_clone.clone()});
    if (mutPlotInfo) showPlotWithInfo(mutPlotInfo, minAF);
}

function write_to_excel(table) {
    var csv = "";
    var cosmRegexp = /id=([0-9]+)/;
    var dbsnpRegexp = /rs=([0-9]+)/;

    table.find("tr").each(function () {
      var sep = "";
      var val = "";
      var db_id = "";
      $(this).find("th").each(function () {
          csv += sep + $(this).text();
          sep = "\t";
      });
      $(this).find("td").each(function () {
          val = $(this).text();
          var links = $(this).find("a");
          for (var a = 0; a < links.length; a++) {
            if (links[a].text.indexOf('COSM') != -1) {
              db_id = links[a].href.match(cosmRegexp);
              val = val.replace('COSM', 'COSM' + db_id[1]);
            }
            if (links[a].text.indexOf('dbSNP') != -1) {
              db_id = links[a].href.match(dbsnpRegexp);
              val = val.replace('dbSNP', 'rs' + db_id[1]);
            }
          }
          csv += sep + val;
          sep = "\t";
      });
      csv += "\n";
    });
    window.URL = window.URL || window.webkiURL;
    var data_type = 'data:application/csv;charset=utf-8,';
    $("#download_mut_table").
    attr("href", data_type + encodeURIComponent(csv)).
    attr("download", "mutations.xls");
}

jQuery(function($) {
    if (!$('#circos_zoom')[0]) return;
    $('#circos_zoom').easyZoom({
        parent: '#circos_plot_div',
        append: false
    });
});

function checkBlacklisted(row, showBlacklisted) {
    for (var c = 0, m = row.cells.length; c < m; c++) {
        var cell = row.cells[c];
        if (cell.attributes.metric && cell.attributes.metric.value == 'VarDict status') {
            if (!showBlacklisted) {  // normal view, hide incidentalome, show else
                if ($(cell).has(".span_status_incidentalome").length > 0)
                    $(row).addClass('row_hidden');
                else
                    $(row).removeClass('row_hidden');
            } else {  // incidentalome view, hide everything else
                if ($(cell).has(".span_status_incidentalome").length == 0)
                    $(row).addClass('row_hidden');
                else
                    $(row).removeClass('row_hidden');
            }
        }
    }
}

function filterMutationsByAF(thresholdValue) {
    minAF = thresholdValue;
    var table_short = $('.table_short#report_table_mutations');
    var table_full = $('.table_full#report_table_mutations');
    var table_actionable = $('#report_table_actionable');
    if (table_short) {
        $(table_short).css('height', '');
        $('.table_short#report_table_mutations tr').each(function() {
            checkAF(this, minAF);
        });
    }
    if (table_full) {
        $(table_full).css('height', '');
        $('.table_full#report_table_mutations tr').each(function() {
            checkAF(this, minAF);
        });
    }
    if (table_actionable) {
        $(table_actionable).css('height', '');
        $('#report_table_actionable tr').each(function() {
            checkAF(this, minAF);
        });
        $('#report_table_actionable tr').each(function() {
            correctRowspan(this);
        });
    }
    if (mutPlotInfo) showPlotWithInfo(mutPlotInfo, minAF);
}

function checkAF(row, minAF) {
      var isKnown = false;
      if (!$('#mut_af_textbox')[0]) return;
      var minActAF = $('#act_min_af')[0].innerText;
      for (var c = 0, m = row.cells.length; c < m; c++) {
          if (row.cells[c].attributes.metric && row.cells[c].attributes.metric.value.indexOf('status') != -1) {
              if (row.cells[c].innerText.indexOf('known') == 0)
                  isKnown = true;
              break;
          }
      }
      for (var c = 0, m = row.cells.length; c < m; c++) {
        if (row.cells[c].attributes.metric && row.cells[c].attributes.number && row.cells[c].attributes.metric.value.indexOf('Freq') != -1) {
            if ((isKnown && row.cells[c].attributes.number.value * 100 >= minActAF) || row.cells[c].attributes.number.value * 100 >= minAF)
                $(row).removeClass('af_less_threshold');
            else $(row).addClass('af_less_threshold');
        }
    }
}
function correctRowspan(row) {
    if (!$(row).find("td:first-child")[0])
        return;
    if ($(row).hasClass("af_less_threshold")) return;

    var rowspan = 1;
    var nextRow = $(row).next("tr");

    while ($(nextRow).find("td:first-child")[0]) {
        if (nextRow.find("td:first-child")[0].attributes.metric.value == "Gene")
            break;
        if (!$(nextRow).hasClass("af_less_threshold"))
            rowspan++;
        nextRow = $(nextRow).next("tr");
    }
    $("td[rowspan]", $(row)).each(function() {
        $(this).attr("rowspan", rowspan);
    });
}

function createSelectVariants() {
    var data = readJsonFromElement($('#mut_parameters_data_json'));
    if (!data) return;

    var selectContainer = $("#select_table_div");
    var backgroundColors = ['#EFF', '#FFE', '#FFEFEE'];
    var tableHead = document.createElement('thead');
    var tableBody = document.createElement('tbody');
    var tableRows = [];
    var maxRow = 0;
    for (var i = 0; i < data.length; i++) {
        maxRow = Math.max(data[i].values.length, maxRow)
    }
    for (var r = 0; r < maxRow; r++) {
        tableRows.push([]);
        for (var i = 0; i < data.length; i++) {
            var td = document.createElement('td');
            tableRows[r].push(td);
        }
    }

    for (var i = 0; i < data.length; i++) {
        values = data[i].values;
        valuesNames = data[i].valuesNames;
        parameter = data[i].parameter;
        parameters.push(parameter);
        var th = document.createElement('th');
        th.innerText = parameter;
        $(tableHead).append(th);

        for (var j = 0; j < values.length; j++) {
            var radioBtn = document.createElement('input');
            var radioBtnId = parameter + j;
            radioBtn.type = "radio";
            radioBtn.name = parameter;
            radioBtn.value = values[j];
            radioBtn.id = radioBtnId;
            if (j == 0)
                $(radioBtn).prop("checked", true);
            if (i == 0 && values[j] != 'all' && values[j] != 'common') {
                var bgColor = backgroundColors[(j - 1) % backgroundColors.length];
                $('#report_table_mutations').find('tr:not(.known) td[metric~="' + values[j].capitalize() + '"]').css("backgroundColor", bgColor);
            }

            var label = document.createElement('label');
            label.htmlFor = radioBtnId;
            label.appendChild(document.createTextNode(valuesNames[j]));
            label.className = "parameter_select";
            $(radioBtn).on("change", function () {
                checkVariantsTable(this.name, this.value);
            });
            td = tableRows[j][i];
            $(td).append(radioBtn);
            $(td).append(label);
        }
    }
    for (var row = 0; row < tableRows.length; row++) {
        var tr = document.createElement('tr');
        for (var cell = 0; cell < tableRows[row].length; cell++) {
            $(tr).append(tableRows[row][cell]);
        }
        $(tableBody).append(tr);
    }
    selectContainer.append(tableHead);
    selectContainer.append(tableBody);
}

function checkVariantsTable(parameter, switchValue) {
    var tables_short = $('.table_short');
    for (var t = 0; t < tables_short.length; t++){
        var table_short = $(tables_short[t]);
        if (table_short[0]) {
            $(table_short).css('height', '');
            $(table_short).find('tbody tr').each(function() {
                checkSamples(this, parameter, switchValue);
            });
        }
    }
}

function checkSamples(row, metric, value) {
    value = value.toLowerCase();
    if (value == 'all') {
        showHideRow(row, metric);
        return;
    }
    if ($(row).attr('class') && $(row).attr('class').indexOf('top_row_tr') != -1)
        return;

    var sampleFound = [];
    var logRatioCols = false;
    var values = value.split(', ');
    for (var c = 0, m = row.cells.length; c < m; c++) {
        var cell = row.cells[c];
        if (cell.attributes.metric) {
            metricName = cell.attributes.metric.value;
            metricValue = cell.innerText;
            if (metricName == metric) {
                if (metricValue.toLowerCase() != value)
                    $(row).addClass(metric + ' unselected_type');
                else showHideRow(row, metric);
            }
            else if (metricName.indexOf('log ratio') != -1) {
                logRatioCols = true;
                if (metricValue) {
                    var isValueSelected = false;
                    for (var v = 0; v < values.length; v++) {
                        if (sampleFound.indexOf(values[v]) == -1 && metricName.toLowerCase().indexOf(values[v]) != -1) {
                            sampleFound.push(values[v]);
                            isValueSelected = true;
                        }
                    }
                    if (!isValueSelected) {
                        $(row).addClass(metric + ' unselected_type');
                        return;
                    }
                }
            }
        }
    }
    if (logRatioCols) {
        if (sampleFound.length != values.length)
            $(row).addClass(metric + ' unselected_type');
        else showHideRow(row, metric);
    }
}

function showHideRow(row, metric) {
    if ($(row).hasClass(metric)) {
        $(row).removeClass(metric);
        for (var p = 0; p < parameters.length; p++)
            if ($(row).hasClass(parameters[p]))
                return;
        $(row).removeClass('unselected_type')
    }
}

function commentMutation(caller) {
    var row = $(caller).parents("tr")[0];
    var gene, transcript, mut, pos, change, effect, significance;
    for (var c = 0; c < row.cells.length; c++) {
        var cell = row.cells[c];
        if (cell.attributes.metric && cell.attributes.metric.value == "Gene") {
            var geneInfo = cell.innerText.split(/\b(\s)/);
            if (geneInfo.length > 1) {
                transcript = geneInfo[0];
                gene = geneInfo[geneInfo.length-1];
            }
            else {
                transcript = "";
                gene = geneInfo[0];
            }
        } else if (cell.attributes.metric && cell.attributes.metric.value == "AA chg") {
            mut = cell.textContent;
        } else if (cell.attributes.metric && cell.attributes.metric.value == "Position") {
            pos = cell.textContent;
        } else if (cell.attributes.metric && cell.attributes.metric.value == "Change") {
            change = cell.textContent;
        } else if (cell.attributes.metric && cell.attributes.metric.value == "Effect") {
            effect = cell.textContent;
        } else if (cell.attributes.metric && cell.attributes.metric.value == "VarDict status") {
            significance = cell.textContent;
        }
    }
    var projectUrl = window.location.href;
    var projectPath = document.getElementById('project_dirpath').innerText;
    document.getElementById('comment_window_text').innerText = 'Leave a comment about mutation ' + mut + ' in ' + gene +
        '. This information is going to be sent to Vlad Saveliev and considered to be added into filtering ' +
        'blacklisting or prioritizing rules.';
    document.getElementById('comment_window_save_btn').onclick=function() {
        var comment = document.getElementById('comment_window_textarea').value;
        if (comment) {
            var data = projectUrl + "," + projectPath + "," + gene + "," + transcript + "," + mut + "," + pos + "," + change +
                "," + effect + "," + significance + "," + comment + "," + getTodayDate() + "\n";
            var php_path = readJsonFromElement($('#comment_php_json'));
            $.post(php_path, {data: data})
             .done(function () {
                  alert('Comment was successfully sent!');
             })
             .error(function () {
                  alert('Error! Comment was not sent.');
             });
        }
        document.getElementById('comment_window').style.display = "none";
    };
    document.getElementById('comment_window_cancel_btn').onclick=function() {
        document.getElementById('comment_window').style.display = "none";
    };

    document.getElementById('comment_window').style.display = "block";

    document.getElementById('comment_window_textarea').focus();

    //document.getElementById("comment_window_textarea").addEventListener("keyup", function(event) {
    //    event.preventDefault();
    //    if (event.keyCode == 13) {
    //        document.getElementById("comment_window_save_btn").click();
    //    }
    //})
}

function getTodayDate(){
    var today = new Date();
    var dd = today.getDate();
    var mm = today.getMonth() + 1;
    var yyyy = today.getFullYear();
    today = dd + '/' + mm + '/' + yyyy;
    return today;
}

function convertHex(hex, opacity){
    hex = hex.replace('#','');
    r = parseInt(hex.substring(0, 2), 16);
    g = parseInt(hex.substring(2, 4), 16);
    b = parseInt(hex.substring(4, 6), 16);

    result = 'rgba('+r+','+g+','+b+','+opacity/100+')';
    return result;
}

String.prototype.capitalize = function() {
    return this.charAt(0).toUpperCase() + this.slice(1);
};

//function extendedClick() {
//    //$('.row_to_hide').toggleClass('row_hidden');
//
//    if (link.html() == 'Full') {
//    } else {
//    }
//}