$(function() {
    var table_low = $('#flagged_low_table_div');
    var table_high = $('#flagged_high_table_div');

    if (msieversion() == 0) {
        table_low.find('table').tableSort();
        table_high.find('table').tableSort();
    }

    if (table_low.find('tr').length > 15 &&
        table_low.find('tr.no_hotspots').length < table_low.find('tr').length - 1) {
        //reduceClick($('#hotspots_low_switch'));
    } else {
        extendClick($('#extend_link'));
    }

    if (table_high.find('tr').length > 15 &&
        table_high.find('tr.no_hotspots').length < table_high.find('tr').length - 1) {
        //reduceClick($('#hotspots_high_switch'));
    } else {
        extendClick($('#extend_link_high'));
    }
});

function extendClick(caller) {
    // Showing full
    var switch_el = $(caller).parent();
    switch_el.html('<a class="dotted-link" id="reduce_link" onclick="reduceClick($(this))">Regions with hotpots & deleterious</a> / <span>all regions</span>');
    var table = switch_el.attr('id').indexOf('low') != -1 ? $('#flagged_low_table_div') : $('#flagged_high_table_div');
    table.find('table').toggleClass('table_short');
}

function reduceClick(caller) {
    // Showing reduced
    var switch_el = $(caller).parent();
    switch_el.html('<span>Regions with hotpots & deleterious</span> / <a class="dotted-link" id="extend_link" onclick="extendClick($(this))">all regions</a>')
    var table = switch_el.attr('id').indexOf('low') != -1 ? $('#flagged_low_table_div') : $('#flagged_high_table_div');
    table.find('table').toggleClass('table_short');
}