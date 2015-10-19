$(function() {
    $("[rel=tooltip]").tooltip({animation: false});

    if (msieversion() == 0) {
        $('table.tableSorter').tableSort();
    }

    if (msieversion() != 0) {
        $(function () {
            //$("td.long_line").style("width: 80px; overflow: clip; text-overflow: clip; white-space: wrap");
            //$("td.long_line span").style("width: 80px; overflow: clip; text-overflow: clip; white-space: wrap");
        });
    }
});