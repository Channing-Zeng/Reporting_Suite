var parameters = [];

$(function() {
    createSelectVariants();
});

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