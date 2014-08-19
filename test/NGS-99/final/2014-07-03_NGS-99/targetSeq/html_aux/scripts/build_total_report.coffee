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

record =
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


DRAGGABLE_COLUMNS = false

RED_HUE = 0
GREEN_HUE = 120
GREEN_HSL = 'hsl(' + GREEN_HUE + ', 80%, 40%)'
CSS_PROP_TO_COLOR = 'background-color'  # color
get_color = (hue) ->
    # lightness = Math.round (Math.pow hue - 75, 2) / 350 + 35
    lightness = 92
    return 'hsl(' + hue + ', 80%, ' + lightness + '%)'


get_meta_tag_contents = (rec) ->
    meta = rec.meta

    if meta? and (a for a of meta).length != 0
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
                    meta_table += "<td>#{toPrettyString(values[db], rec.metric.unit)}</td>"

                meta_table += '</tr>\n'
            meta_table += '</table>\n'

            return "class=\"meta_info_span tooltip-meta\" rel=\"tooltip\" title=\"#{meta_table}\""
    else
        return "class=\"meta_info_span tooltip-meta\" rel=\"tooltip\""


get_metric_name_html = (rec) ->
    if rec.metric.short_name
        metricName = rec.metric.short_name
    else
        metricName = rec.metric.name
    if rec.metric.description
        return "<a class=\"tooltip-link\" rel=\"tooltip\" title=\"#{rec.metric.description}\">
            #{metricName}
        </a>"
    else
        return metricName


calc_records_cell_contents = (records, font) ->
    for rec in records
        value = rec.value
        num_html = ''

        if not value? or value == ''
            rec.cell_contents = '-'

        else
            if typeof value == 'number'
                rec.num = value
                rec.cell_contents = toPrettyString value, rec.metric.unit
                num_html = toPrettyString value

            else if /^-?.?[0-9]/.test(value)
                result = /([0-9\.]+)(.*)/.exec value
                rec.num = parseFloat result[1]
                rec.cell_contents = toPrettyString(rec.num, rec.metric.unit) + result[2]
                num_html = toPrettyString(rec.num)
            else
                rec.cell_contents = value

        # Max frac width of column
        rec.frac_width = $.fn.intPartTextWidth num_html, font


calc_cell_contents = (report, font) ->
    max_frac_widths_by_metric = {}
    min_val_by_metric = {}
    max_val_by_metric = {}

    for sampleReport in report.sample_reports
        calc_records_cell_contents sampleReport.records, font
        for rec in sampleReport.records
            if not (rec.metric.name of max_frac_widths_by_metric)
                max_frac_widths_by_metric[rec.metric.name] = rec.frac_width
            else if rec.frac_width > max_frac_widths_by_metric[rec.metric.name]
                max_frac_widths_by_metric[rec.metric.name] = rec.frac_width

            # Max and min value (for heatmap)
            if rec.num?
                if not (rec.metric.name of min_val_by_metric)
                    min_val_by_metric[rec.metric.name] = rec.num
                else if min_val_by_metric[rec.metric.name] > rec.num
                    min_val_by_metric[rec.metric.name] = rec.num

                if not (rec.metric.name of max_val_by_metric)
                    max_val_by_metric[rec.metric.name] = rec.num
                else if max_val_by_metric[rec.metric.name] < rec.num
                    max_val_by_metric[rec.metric.name] = rec.num

    for sampleReport in report.sample_reports
        for rec in sampleReport.records
            # Padding based on frac width
            if rec.frac_width?
                rec.right_shift = max_frac_widths_by_metric[rec.metric.name] - rec.frac_width

                if rec.right_shift != 0
                    a = 0

            # Color heatmap
            if rec.num?
                max = max_val_by_metric[rec.metric.name]
                min = min_val_by_metric[rec.metric.name]

                maxHue = GREEN_HUE
                minHue = RED_HUE
                if rec.metric.quality == 'Less is better'
                    maxHue = RED_HUE
                    minHue = GREEN_HUE

                if max == min
#                    rec.color = get_color GREEN_HUE
                else
                    k = (maxHue - minHue) / (max - min)
                    hue = Math.round minHue + (rec.num - min) * k
                    rec.color = get_color hue


reporting.buildCommonRecords = (common_records) ->
    if common_records
        calc_records_cell_contents common_records, $('#report').css 'font'

        table = "<table cellspacing=\"0\" class=\"common_table\" id=\"common_table\">"
        for rec in common_records
            table += "\n<tr><td>
                    <span class='metric_name'>#{get_metric_name_html(rec)}:</span>
                    #{rec.cell_contents}
                  </td></tr>"
        table += "\n</table>\n"

        $('#report').append table


reporting.buildTotalReport = (report, columnOrder) ->
    if report.name?
        $('#report').append "<h3 class='table_name' style='margin: 0px 0 5px 0'>#{report.name}</h3>"

    calc_cell_contents report, $('#report').css 'font'

    table = "<table cellspacing=\"0\"
                    class=\"report_table #{if DRAGGABLE_COLUMNS then 'draggable' else ''} fix-align-char\"
                    id=\"report_table_#{report.name}\">"
    table += "\n<tr class=\"top_row_tr\">"
    table += "<td class=\"top_left_td left_column_td\">
                    <span>Sample</span>
              </td>"

    for recNum in [0...report.sample_reports[0].records.length]
        pos = columnOrder[recNum]
        rec = report.sample_reports[0].records[pos]
        table += "<td class='second_through_last_col_headers_td' position='#{pos}'>
             #{if DRAGGABLE_COLUMNS then '<span class=\'drag_handle\'><span class=\'drag_image\'></span></span>' else ''}
             <span class='metricName'>#{get_metric_name_html(rec)}</span>
        </td>"

    for sampleReport in report.sample_reports
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

            table += "<td metric=\"#{rec.metric.name}\"
                          style=\"#{CSS_PROP_TO_COLOR}: #{rec.color}\"
                          class='number'
                          quality=\"#{rec.metric.quality}\""
            if rec.num? then table += ' number="' + rec.value + '">'
            if rec.right_shift?
                padding = "margin-left: #{rec.right_shift}px; margin-right: -#{rec.right_shift}px;"
            else
                padding = ""
            table += "<a style=\"#{padding}\"
                          #{get_meta_tag_contents(rec)}>#{rec.cell_contents}
                      </a>
                      </td>"
        table += "</tr>"
    table += "\n</table>\n"

    $('#report').append table


set_legend = ->
    legend = '<span>'
    step = 6
    for hue in [RED_HUE..GREEN_HUE] by step

        legend += "<span style=\"#{CSS_PROP_TO_COLOR}: #{get_color hue}\">"

        switch hue
            when RED_HUE              then legend += 'w'
            when RED_HUE   +     step then legend += 'o'
            when RED_HUE   + 2 * step then legend += 'r'
            when RED_HUE   + 3 * step then legend += 's'
            when RED_HUE   + 4 * step then legend += 't'
            when GREEN_HUE - 3 * step then legend += 'b'
            when GREEN_HUE - 2 * step then legend += 'e'
            when GREEN_HUE -     step then legend += 's'
            when GREEN_HUE            then legend += 't'
            else                           legend += '.'
        legend += "</span>"
    legend += "</span>"
    $('#report_legend').append legend


$.fn._splitDot_partTextWidth = (html, font, part_type) ->  # part_type = 'int'|'frac'
    parts = html.split '.'

    if part_type == 'frac'
        if parts.length < 2
            return 0
        else
            frac_part = '.' + parts[1]

    else if part_type == 'int'
        frac_part = parts[0]

    if (!$.fn.fracPartTextWidth.fakeEl)
        $.fn.fracPartTextWidth.fakeEl = $('<span>').hide().appendTo document.body

    $.fn.fracPartTextWidth.fakeEl.html frac_part
    $.fn.fracPartTextWidth.fakeEl.css 'font', font
    return $.fn.fracPartTextWidth.fakeEl.width()


$.fn.fracPartTextWidth = (html, font) ->
    $.fn._splitDot_partTextWidth html, font, 'frac'


$.fn.intPartTextWidth = (html, font) ->
    $.fn._splitDot_partTextWidth html, font, 'int'


$.fn.textWidth = (text, font) ->
    if (!$.fn.textWidth.fakeEl)
        $.fn.textWidth.fakeEl = $('<span>').hide().appendTo document.body

    $.fn.textWidth.fakeEl.html text
    $.fn.textWidth.fakeEl.css 'font', font
    return $.fn.textWidth.fakeEl.width()


String.prototype.trunc = (n) ->
    this.substr(0, n - 1) + (this.length > n ? '&hellip;': '')



#postprocess_cells = ->
#    processes_metrics = []
#
#    $(".report_table td[number]").each ->
#        metricName = $(this).attr 'metricName'
#
#        if !(metricName in processes_metrics)
#            processes_metrics.push metricName
#            console.log metricName
#
#            quality = $(this).attr 'quality'
#            all_cells = $('.report_table').find "td[metricName=\"#{metricName}\"][number]"
#            all_numbers = ($(cell).attr 'number' for cell in all_cells)
#
#            set_heatmap all_cells, all_numbers, quality
#
#            set_offset all_cells, all_numbers, metricName