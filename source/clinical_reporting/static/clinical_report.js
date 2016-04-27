$(function() {
    $("[rel=tooltip]").tooltip({animation: false});

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

            if (table_short.find('tr').length > 1 &&
                table_full.find('tr').length > 15) {
                reduceClick('extend_link_' + switch_id);
            } else {
                extendClick('extend_link_' + switch_id);
            }
        }
    }
    if ($('#mut_af_slider')[0]) filterMutationsByAF($('#mut_af_slider')[0].value);

    $('#variants_table_controls').width($('#report_table_mutations').width() - 5);
    $('#download_mut_table').show();
    //if (msieversion() == 0) {
    //    $('table.tableSorter.table_short').tableSort();
    //}
});

function extendClick(switch_id) {
    if (switch_id[0].id) switch_id = switch_id[0].id;
    var showSilent = switch_id.search('silent') != -1;
    // Showing full
    switch_id = switch_id.split("_");
    var table_id = switch_id[switch_id.length - 1];
    var switch_el = $('#' + table_id + '_switch');
    if (table_id == 'seq2c') {
      switch_el.html('<a class="dotted-link" id="reduce_link_' + table_id + '" onclick="reduceClick($(this))">' + key_or_target + ' genes</a> / <span>all genes</span>')
    }
    else {
        switchElContent = '<a class="dotted-link" id="reduce_link_' + table_id + '" onclick="reduceClick($(this))">known, likely</a> / ';
        if (showSilent) {
            switchElContent += '<a class="dotted-link" id="extend_link_' + table_id + '" onclick="extendClick($(this))">+ unknown</a> / ';
            switchElContent += '<span>+ silent</span>';
        }
        else {
            switchElContent += '<span>+ unknown</span> / ';
            switchElContent += '<a class="dotted-link" id="extend_link_silent_' + table_id + '" onclick="extendClick($(this))">+ silent</a>';
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
                checkSilent(this, showSilent);
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
}

function reduceClick(switch_id) {
  if (switch_id[0].id) switch_id = switch_id[0].id;
    switch_id = switch_id.split("_");
    var table_id = switch_id[switch_id.length - 1];
    // Showing reduced
    var switch_el = $('#' + table_id + '_switch');
    if (table_id == 'seq2c') {
      switch_el.html('<span>'  + key_or_target + ' genes</span> / <a class="dotted-link" id="extend_link_' + table_id + '" onclick="extendClick($(this))">all genes</a>')
    }
    else {
      switch_el.html('<span>known, likely</span> / <a class="dotted-link" id="extend_link_' + table_id + '" onclick="extendClick($(this))">+ unknown</a> ' +
          '/ <a class="dotted-link" id="extend_link_silent_' + table_id + '" onclick="extendClick($(this))">+ silent</a>')
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
    $('#circos_zoom').easyZoom({
        parent: '#circos_plot_div',
        append: false
    });
});

function checkSilent(row, showSilent) {
    for (var c = 0, m = row.cells.length; c < m; c++) {
        var cell = row.cells[c];
        if (cell.attributes.metric && cell.attributes.metric.value == 'VarDict status') {
            if (!showSilent && cell.innerText.search('silent') != -1)
                $(row).addClass('row_hidden');
            else $(row).removeClass('row_hidden');
        }
    }
}

function filterMutationsByAF(thresholdValue) {
    minAF = thresholdValue;
    var table_short = $('.table_short#report_table_mutations');
    var table_full = $('.table_full#report_table_mutations');
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
}

function checkAF(row, minAF) {
      for (var c = 0, m = row.cells.length; c < m; c++) {
        if (row.cells[c].attributes.metric && row.cells[c].attributes.number && row.cells[c].attributes.metric.value.indexOf('Freq') != -1) {
            if (row.cells[c].attributes.number.value * 100 < minAF)
                $(row).addClass('af_less_threshold');
            else $(row).removeClass('af_less_threshold');
        }
    }
}
function showVariantsByType(switch_id) {
    if (switch_id[0].id) switch_id = switch_id[0].id;
    // Showing full
    switch_id = switch_id.split("_");
    var switchValue = switch_id[switch_id.length-1];
    var switch_el = $('#variants_group_switch_type');
    var parameter = 'Type';
    if (switchValue == 'all')
        switchElContent = '<span>all mutations</span> / ';
    else
        switchElContent = '<a class="dotted-link" id="show_variants_all" onclick="showVariantsByType($(this))"> all mutations</a> / ';
    if (switchValue == 'plasma')
        switchElContent += '<span> only plasma samples</span> / ';
    else
        switchElContent += '<a class="dotted-link" id="show_variants_plasma" onclick="showVariantsByType($(this))"> only plasma samples</a> / ';
    if (switchValue == 'tissue')
        switchElContent += '<span> only tissue samples</span> / ';
    else
        switchElContent += '<a class="dotted-link" id="show_variants_tissue" onclick="showVariantsByType($(this))"> only tissue samples</a> / ';

    switch_el.html(switchElContent);
    checkVariantsTable(parameter, switchValue);
}

function showVariantsBySensitivity(switch_id) {
    if (switch_id[0].id) switch_id = switch_id[0].id;
    // Showing full
    switch_id = switch_id.split("_");
    var switchValue = switch_id[switch_id.length-1];
    var switch_el = $('#variants_group_switch_sens');
    var parameter = 'Sensitivity';
    if (switchValue == 'all')
        switchElContent = '<span>all mutations</span> / ';
    else
        switchElContent = '<a class="dotted-link" id="show_variants_all" onclick="showVariantsBySensitivity($(this))"> all mutations</a> / ';
    if (switchValue == 'sensitive')
        switchElContent += '<span> only sensitive samples</span> / ';
    else
        switchElContent += '<a class="dotted-link" id="show_variants_sensitive" onclick="showVariantsBySensitivity($(this))"> only sensitive samples</a> / ';
    if (switchValue == 'resistant')
        switchElContent += '<span> only resistant samples</span> / ';
    else
        switchElContent += '<a class="dotted-link" id="show_variants_resistant" onclick="showVariantsBySensitivity($(this))"> only resistant samples</a> / ';
    if (switchValue == 'common')
        switchElContent += '<span> only common mutations</span> / ';
    else
        switchElContent += '<a class="dotted-link" id="show_variants_common" onclick="showVariantsBySensitivity($(this))"> only common mutations</a> / ';
    if (switchValue == 'common') switchValue = 'resistant, sensitive';

    switch_el.html(switchElContent);
    checkVariantsTable(parameter, switchValue);
}

function checkVariantsTable(parameter, switchValue) {
    var table_short = $('.table_short#report_table_mutations');
    var table_full = $('.table_full#report_table_mutations');
    if (table_full[0]) {
        $(table_full).css('height', '');
        $('.table_full#report_table_mutations tr').each(function() {
            checkSamples(this, parameter, switchValue);
        });
    }
    else {
        $(table_short).css('height', '');
        $('.table_short#report_table_mutations tr').each(function() {
            checkSamples(this, parameter, switchValue);
        });
    }
}

function checkSamples(row, metric, value) {
    for (var c = 0, m = row.cells.length; c < m; c++) {
        var cell = row.cells[c];
        if (value == 'all') $(row).removeClass('unselected_type');
        else {
            if (cell.attributes.metric && cell.attributes.metric.value == metric) {
                if (cell.innerText.toLowerCase() != value)
                    $(row).addClass('unselected_type');
                else $(row).removeClass('unselected_type');
            }
        }
    }
}
//function extendedClick() {
//    //$('.row_to_hide').toggleClass('row_hidden');
//
//    if (link.html() == 'Full') {
//    } else {
//    }
//}