String.prototype.trunc = (n) ->
    this.substr(0, n - 1) + (this.length > n ? '&hellip;': '')

report =
    name: ''
    fpath: ''
    metrics: []

metric =
    name: ''
    quality: ''
    value: ''
    meta: ''


get_meta_span = (meta) ->
    if meta? and (k for own k of meta).length isnt 0
        meta_table = '<table class=\'qc_meta_table\'>\n<tr><td></td>'
        for novelty, values of meta
            meta_table += '<td>' + novelty + '</td>'
        meta_table += '</tr>\n'

        for novelty, values of meta
            dbs = (db for db, val of values)
            break

        for db in dbs
            meta_table += '<tr>' + '<td>' + db + '</td>'
            for novelty, values of meta
                meta_table += '<td>' + values[db] + '</td>'
            meta_table += '</tr>\n'
        meta_table += '</table>\n'

        return 'class="meta_info_span" rel="tooltip" title="' + meta_table + '"'
    else
        return ''


RED_HUE = 0
GREEN_HUE = 120
GREEN_HSL = 'hsl(' + GREEN_HUE + ', 80%, 40%)'


reporting.buildTotalReport = (report, columnOrder, date, glossary) ->
    $('#report_date').html('<p>' + date + '</p>');

    table = '<table cellspacing="0" class="report_table draggable" id="main_report_table">'
    table += '<tr class="top_row_tr">
                      <td id="top_left_td" class="left_column_td">
                        <span>Sample</span>
                      </td>'

    for metricNum in [0...report[0].metrics.length]
        pos = columnOrder[metricNum]
        metric = report[0].metrics[pos]
        table +=
            '<td class="second_through_last_col_headers_td" position="' + pos + '">' +
            '<span class="drag_handle"><span class="drag_image"></span></span>' +
            '<span class="metric_name">' +
            nbsp(addTooltipIfDefinitionExists(glossary, metric.name), metric.name) +
            '</span>
                               </td>'

    for sampleReport in report
        sampleName = sampleReport.name
        if sampleReport.name.length > 30
            sampleName = '<span class="tooltip-link" rel="tooltip" title="' + sampleName + '">' +
            sampleName.trunc(80) + '</span>'

        table += '<tr><td class="left_column_td"><a class="sample_name" href="html_aux/' + sampleName + '.html">' + sampleName + '</a></td>'

        for metricNum in [0...sampleReport.metrics.length]
            pos = columnOrder[metricNum]
            metric = sampleReport.metrics[pos]
            quality = metric.quality
            value = metric.value
            meta_span = get_meta_span metric.meta

            table += '<td metric="' + metric.name + '" quality="' + metric.quality + '"'

            if not value? or value == ''
                cell_contents = '-'

            else
                if typeof value == 'number'
                    num = value
                    cell_contents = toPrettyString(value)
                else
                    result = /([0-9\.]+)(.*)/.exec value
                    num = parseFloat result[1]
                    cell_contents = toPrettyString(num) + result[2]

                if num?
                    table += ' number="' + value + '">'

            table += '<span ' + meta_span + '>' + cell_contents + '</span></td>'
        table += '</tr>'
    table += '</table>'

    $('#report').append table
    $('#report_legend').append legend
    add_heatmap()

    legend = '<span>';
    step = 6;
    for hue in [RED_HUE..GREEN_HUE] by step
        lightness = (Math.pow(hue - 75, 2)) / 350 + 35
        legend += '<span style="color: hsl(' + hue + ', 80%, ' + lightness + '%);">'

        switch hue
            when RED_HUE              then legend += 'w'
            when RED_HUE + step       then legend += 'o'
            when RED_HUE + 2 * step   then legend += 'r'
            when RED_HUE + 3 * step   then legend += 's'
            when RED_HUE + 4 * step   then legend += 't'
            when GREEN_HUE - 3 * step then legend += 'b'
            when GREEN_HUE - 2 * step then legend += 'e'
            when GREEN_HUE - step     then legend += 's'
            when GREEN_HUE            then legend += 't'
            else
                legend += '.'
        legend += '</span>'
    legend += '</span>'


add_heatmap = ->
    $(".report_table td[number]").each ->
        metric = $(this).attr 'metric'
        quality = $(this).attr 'quality'
        other_cells = $('.report_table').find 'td[metric="' + metric + '"][number]'
        orther_numbers = ($(cell).attr 'number' for cell in other_cells)

        min = Math.min.apply(null, orther_numbers)
        max = Math.max.apply(null, orther_numbers)

        maxHue = GREEN_HUE
        minHue = RED_HUE

        if quality == 'Less is better'
            maxHue = RED_HUE
            minHue = GREEN_HUE

        if max == min
            $(other_cells).css('color', GREEN_HSL)
        else
            k = (maxHue - minHue) / (max - min)
            hue = 0
            lightness = 0
            other_cells.each (i) ->
                number = orther_numbers[i]
                hue = Math.round minHue + (number - min) * k
                lightness = Math.round (Math.pow(hue - 75, 2)) / 350 + 35
                # $(this).css('color', 'hsl(' + hue + ', 80%, 35%)');
                $(this).css('color', 'hsl(' + hue + ', 80%, ' + lightness + '%)')

                $('#report_legend').show('fast') if orther_numbers.length > 1
