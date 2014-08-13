String.prototype.trunc = (n) ->
    this.substr(0, n - 1) + (this.length > n ? '&hellip;': '')

report =
    sample:
        name: ''
        phenotype: ''
        bam: ''
        bed: ''
        vcf_by_caller:
            name: ''
            summary_qc_rep_fpaths: []
            anno_vcf_fpaths: {}
            anno_filt_vcf_fpaths: {}
    fpath: ''
    link: ''
    records: []

records =
    metric: null
    value: ''
    meta: null

metric =
    name: ''
    short_name: ''
    description: ''
    quality: ''
    presision: 0
    type: null


get_meta_tag_contents = (rec) ->
    metric = rec.metric
    meta = rec.meta

    if meta? and meta.length isnt 0
        return "class=\"meta_info_span tooltip-meta\" rel=\"tooltip\""

        if typeof meta is 'string'
            return "class=\"meta_info_span tooltip-meta\" rel=\"tooltip\" title=\"#{meta}\""

        else  # qc
            (k for own k of meta).length isnt 0
            meta_table = '<table class=\'qc_meta_table\'>\n<tr><td></td>'
            for novelty, values of meta when novelty isnt 'all'
                meta_table += "<td>#{novelty}</td>"
    #        for novelty, values of meta when novelty is 'all'
    #            meta_table += "<td>#{novelty}</td>"
            meta_table += '</tr>\n'

            for novelty, values of meta
                dbs = (db for db, val of values when db isnt 'average')
                dbs.push 'average'
                break

            for db in dbs
                meta_table += "<tr><td>#{db}</td>"
                for novelty, values of meta when novelty isnt 'all'
                    meta_table += "<td>#{toPrettyString(values[db], metric.unit)}</td>"
    #            for novelty, values of meta when novelty is 'all'
    #                meta_table += "<td>#{toPrettyString(values[db], metric.unit)}</td>"

                meta_table += '</tr>\n'
            meta_table += '</table>\n'

            return "class=\"meta_info_span tooltip-meta\" rel=\"tooltip\" title=\"#{meta_table}\""


RED_HUE = 0
GREEN_HUE = 120
GREEN_HSL = 'hsl(' + GREEN_HUE + ', 80%, 40%)'


reporting.buildTotalReport = (report, columnOrder, date) ->
    $('#report_date').html('<p>' + date + '</p>');

    table = "<table cellspacing=\"0\" class=\"report_table draggable fix-align-char\" id=\"main_report_table\">"
    table += "\n<tr class=\"top_row_tr\">"
    table += "<td id=\"top_left_td\" class=\"left_column_td\">
                <span>Sample</span>
              </td>"

    for recNum in [0...report[0].records.length]
        pos = columnOrder[recNum]
        rec = report[0].records[pos]
        metric = rec.metric
        if metric.description
            metric_html = "<a class=\"tooltip-link\" rel=\"tooltip\" title=\"#{metric.description}\">
                #{metric.short_name}
            </a>"
        else
            if metric.short_name is undefined
                metric_html = metric.name
            else
                metric_html = metric.short_name

        table += "<td class='second_through_last_col_headers_td' position='#{pos}'>
             <span class='drag_handle'><span class='drag_image'></span></span>
             <span class='metric_name'>#{metric_html}</span>
        </td>"

    for sampleReport in report
        sampleName = sampleReport.sample.name
        sampleLink = sampleReport.link
        if sampleName.length > 30
            sampleName = "<span title=\"#{sampleName}\">#{sampleName.trunc(80)}</span>"

        table += "\n<tr>
            <td class=\"left_column_td\">
                <a class=\"sample_name\" href=\"#{sampleLink}\">#{sampleName}</a>
            </td>"

        for recNum in [0...sampleReport.records.length]
            pos = columnOrder[recNum]
            rec = sampleReport.records[pos]
            metric = rec.metric
            value = rec.value

            if not value? or value == ''
                rec.cell_contents = '-'

            else
                if typeof value == 'number'
                    rec.num = value
                    rec.cell_contents = toPrettyString(value, metric.unit)

                else if /^-?.?[0-9]/.test(value)
                    result = /([0-9\.]+)(.*)/.exec value
                    rec.num = parseFloat result[1]
                    rec.cell_contents = toPrettyString(rec.num, metric.unit) + result[2]

                else
                    rec.cell_contents = value

            table += "<td metric=\"#{metric.name}\" quality=\"#{metric.quality}\""

            if rec.num?
                table += ' number="' + rec.value + '">'

            table += "<a #{get_meta_tag_contents(rec)}>#{rec.cell_contents}</a></td>"
        table += "</tr>"
    table += "\n</table>\n"

    $('#report').append table
    $('#report_legend').append legend
#    $('#report_table').alignColumn [1...report[0].records.length]
    add_heatmap()

    legend = '<span>'
    step = 6
    for hue in [RED_HUE..GREEN_HUE] by step
        lightness = (Math.pow(hue - 75, 2)) / 350 + 35
        legend += "<span style=\"color: hsl(#{hue}, 80%, #{lightness}%);\">"

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
            else                           legend += '.'
        legend += "</span>"
    legend += "</span>"


css_property_to_color = 'background-color'  # color


get_color = (hue) ->
#    lightness = Math.round (Math.pow hue - 75, 2) / 350 + 75
#    lightness = Math.round (Math.pow hue - 75, 2) / 350 + 35
    lightness = 85
    return 'hsl(' + hue + ', 80%, ' + lightness + '%)'


add_heatmap = ->
    processes_metrics = []

    $(".report_table td[number]").each ->
        metric = $(this).attr 'metric'
        if !(metric in processes_metrics)
            quality = $(this).attr 'quality'
            other_cells = $('.report_table').find "td[metric=\"#{metric}\"][number]"
            orther_numbers = ($(cell).attr 'number' for cell in other_cells)

            min = Math.min.apply null, orther_numbers
            max = Math.max.apply null, orther_numbers

            maxHue = GREEN_HUE
            minHue = RED_HUE

            if quality == 'Less is better'
                maxHue = RED_HUE
                minHue = GREEN_HUE

    #        if (num for num in orther_numbers when num < 10).length > 0
    #            $(other_cells).addClass 'float'
    #        else
    #            $(other_cells).addClass 'integer'

            if max == min
                $(other_cells).css css_property_to_color, get_color GREEN_HUE
            else
                k = (maxHue - minHue) / (max - min)
                other_cells.each (i) ->
                    number = orther_numbers[i]
                    hue = Math.round minHue + (number - min) * k
                    $(this).css css_property_to_color, get_color hue

                    $('#report_legend').show 'fast' if orther_numbers.length > 1

            # Decimal point alignment
            left_part_width = 0
            for cell in $(other_cells)
                parts = $(cell).text().split '.'
                if (parts.length > 1)
                    left_part =  '.' + parts[-1]
                else
                    left_part = ''
                width = $.fn.textWidth(left_part)

                left_part_width = width if width > left_part_width

            if left_part_width > 0
                Console.log left_part_width

            for cell in other_cells
                offset = $(cell).offset()
                offset.left += left_part_width
                $(cell).offset(offset)

            processes_metrics.push metric


$.fn.textWidth = (text) ->
    if (!$.fn.textWidth.fakeEl)
        $.fn.textWidth.fakeEl = $('<span>').hide().appendTo(document.body)

    $.fn.textWidth.fakeEl.text(text)

    $.fn.textWidth.fakeEl.css('font', this.css('font'))

    return $.fn.textWidth.fakeEl.width()