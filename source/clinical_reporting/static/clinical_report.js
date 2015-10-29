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

var table_short_clone = null;
var table_full_clone = null;

$(function() {
    var table_short = $('.table_short');
    var table_full = $('.table_full');
    table_short_clone = table_short.clone();
    table_full_clone = table_full.clone();

    table_short.remove();
    table_full.remove();

    if (table_short.find('tr').length > 1 &&
        table_full.find('tr').length > 15) {
        reduceClick();
    } else {
        extendClick();
    }

    //if (msieversion() == 0) {
    //    $('table.tableSorter.table_short').tableSort();
    //}
});

function extendClick() {
    // Showing full
    var switch_el = $('#variants_switch');
    switch_el.html('<a class="dotted-link" id="reduce_link" onclick="reduceClick($(this))">Likely pathogenic</a> / <span>all variants</span>')
    var var_table_div = $('#variants_table_div');
    $('.table_short').remove();
    var_table_div.prepend(table_full_clone);
    $('.table_full').show();
    if (msieversion() == 0) {
        $('table.tableSorter.table_full').tableSort();
    }
    table_full_clone = table_full_clone.clone();
}

function reduceClick() {
    // Showing reduced
    var switch_el = $('#variants_switch');
    switch_el.html('<span>Likely pathogenic</span> / <a class="dotted-link" id="extend_link" onclick="extendClick($(this))">all variants</a>')
    var var_table_div = $('#variants_table_div');
    $('.table_full').remove();
    var_table_div.prepend(table_short_clone);
    $('.table_short').show();
    if (msieversion() == 0) {
        $('table.tableSorter.table_short').tableSort();
    }
    table_short_clone = table_short_clone.clone();
}

//function extendedClick() {
//    //$('.row_to_hide').toggleClass('row_hidden');
//
//    if (link.html() == 'Full') {
//    } else {
//    }
//}