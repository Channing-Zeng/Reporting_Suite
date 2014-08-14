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

metric_name =
    name: ''
    short_name: ''
    description: ''
    quality: ''
    presision: 0
    type: null


get_meta_tag_contents = (rec) ->
    metric_name = rec.metric
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
            meta_table += '</tr>\n'

            for novelty, values of meta
                dbs = (db for db, val of values when db isnt 'average')
                dbs.push 'average'
                break

            for db in dbs
                meta_table += "<tr><td>#{db}</td>"
                for novelty, values of meta when novelty isnt 'all'
                    meta_table += "<td>#{toPrettyString(values[db], metric_name.unit)}</td>"

                meta_table += '</tr>\n'
            meta_table += '</table>\n'

            return "class=\"meta_info_span tooltip-meta\" rel=\"tooltip\" title=\"#{meta_table}\""


RED_HUE = 0
GREEN_HUE = 120
GREEN_HSL = 'hsl(' + GREEN_HUE + ', 80%, 40%)'

DRAGGABLE_COLUMNS = false


reporting.buildTotalReport = (report, columnOrder, date) ->
    $('#report_date').html('<p>' + date + '</p>');

    table = "<table cellspacing=\"0\" class=\"report_table #{if DRAGGABLE_COLUMNS then 'draggable' else ''} fix-align-char\" id=\"main_report_table\">"
    table += "\n<tr class=\"top_row_tr\">"
    table += "<td id=\"top_left_td\" class=\"left_column_td\">
                <span>Sample</span>
              </td>"

    for recNum in [0...report[0].records.length]
        pos = columnOrder[recNum]
        rec = report[0].records[pos]
        metric_name = rec.metric
        if metric_name.description
            metric_html = "<a class=\"tooltip-link\" rel=\"tooltip\" title=\"#{metric_name.description}\">
                #{metric_name.short_name}
            </a>"
        else
            if metric_name.short_name is undefined
                metric_html = metric_name.name
            else
                metric_html = metric_name.short_name

        table += "<td class='second_through_last_col_headers_td' position='#{pos}'>
             #{if DRAGGABLE_COLUMNS then '<span class=\'drag_handle\'><span class=\'drag_image\'></span></span>' else ''}
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
            metric_name = rec.metric
            value = rec.value

            if not value? or value == ''
                rec.cell_contents = '-'

            else
                if typeof value == 'number'
                    rec.num = value
                    rec.cell_contents = toPrettyString(value, metric_name.unit)

                else if /^-?.?[0-9]/.test(value)
                    result = /([0-9\.]+)(.*)/.exec value
                    rec.num = parseFloat result[1]
                    rec.cell_contents = toPrettyString(rec.num, metric_name.unit) + result[2]

                else
                    rec.cell_contents = value

            table += "<td metric_name=\"#{metric_name.name}\" class='number' quality=\"#{metric_name.quality}\""

            if rec.num?
                table += ' number="' + rec.value + '">'

            table += "<a #{get_meta_tag_contents(rec)}>#{rec.cell_contents}</a></td>"
        table += "</tr>"
    table += "\n</table>\n"

    $('#report').append table

    postprocess_cells()


postprocess_cells = ->
    CSS_PROP_TO_COLOR = 'background-color'  # color

    get_color = (hue) ->
        # lightness = Math.round (Math.pow hue - 75, 2) / 350 + 75
        # lightness = Math.round (Math.pow hue - 75, 2) / 350 + 35
        lightness = 85
        return 'hsl(' + hue + ', 80%, ' + lightness + '%)'

    set_heatmap = (all_cells, all_numbers, quality) ->
        min = Math.min.apply null, all_numbers
        max = Math.max.apply null, all_numbers

        maxHue = GREEN_HUE
        minHue = RED_HUE

        if quality == 'Less is better'
            maxHue = RED_HUE
            minHue = GREEN_HUE

        if max == min
            $(all_cells).css CSS_PROP_TO_COLOR, get_color GREEN_HUE
        else
            k = (maxHue - minHue) / (max - min)
            all_cells.each (i) ->
                number = all_numbers[i]
                hue = Math.round minHue + (number - min) * k
                $(this).css CSS_PROP_TO_COLOR, get_color hue

                $('#report_legend').show 'fast' if all_numbers.length > 1

    set_offset = (all_cells, all_numbers, metric_name) ->
        if (num for num in all_numbers when isFractional num).length == 0
            console.log 'Integer: ' + metric_name
        else
            console.log 'Decimal: ' + metric_name

            # Floating point: decimal point alignment
            left_part_widths = []
            max_left_part_width = 0
            for cell in $(all_cells)
                text = $(cell).text()
                parts = text.split '.'
                left = parts[parts.length - 1]

                left_part = if left? and left != 'undefined' then '.' + left else left_part = ''

                width = $.fn.textWidth left_part, $(cell).css 'font'
                left_part_widths.push width
                max_left_part_width = width if width > max_left_part_width

                # if (w for w in left_part_widths when w != max_left_part_width).length > 0
                #    console.log 'max_left_part_width: ' + max_left_part_width

            for i in [0...left_part_widths.length]
                left_part_width = left_part_widths[i]
                if max_left_part_width != left_part_width
                    cell = $(all_cells[i])
                    a = cell.children('a')
                    offset = a.offset()
                    offset.left -= max_left_part_width - left_part_width
                    a.offset offset

    set_legend = ->
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
        $('#report_legend').append legend

    $.fn.textWidth = (text, font) ->
        if (!$.fn.textWidth.fakeEl)
            $.fn.textWidth.fakeEl = $('<span>').hide().appendTo document.body

        $.fn.textWidth.fakeEl.text text
        $.fn.textWidth.fakeEl.css 'font', font

        return $.fn.textWidth.fakeEl.width()

    processes_metrics = []

    $(".report_table td[number]").each ->
        metric_name = $(this).attr 'metric_name'

        if !(metric_name in processes_metrics)
            processes_metrics.push metric_name
            console.log metric_name

            quality = $(this).attr 'quality'
            all_cells = $('.report_table').find "td[metric_name=\"#{metric_name}\"][number]"
            all_numbers = ($(cell).attr 'number' for cell in all_cells)

            set_heatmap all_cells, all_numbers, quality

            set_offset all_cells, all_numbers, metric_name

    set_legend()
