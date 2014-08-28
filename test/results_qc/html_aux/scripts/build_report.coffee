showPlotWithInfo = (info) ->
    newSeries = []
    newColors = []

    $('#legend-placeholder').find 'input:checked'.each ->
        number = $(this).attr 'name'
        if number and info.series && info.series.length > 0
            for i in [i...info.series.length]
                series = info.series[i]
                break if series.number != number

            if i <= info.series.length
                newSeries.push(series)
                newColors.push(series.color)
            else
                console.log('no series with number ' + number)

        if newSeries.length == 0
            newSeries.push
                data: []

            newColors.push '#FFF'

        info.showWithData(newSeries, newColors)


recoverOrderFromCookies = (report_name) ->
    return null unless navigator.cookieEnabled

    orderString = readCookie report_name + '_order'
    return null unless orderString

    columnOrder = []
    fail = false

    for val in orderString.split(' ')
        val = parseInt val
        if isNaN val
            fail = true
        else
            columnOrder.push val

    return null if fail

    return columnOrder


readJson = (what) ->
    result
    try
        result = JSON.parse $('#' + what + '-json').html()
    catch e
        result = null

    return result


totalReportData =
    date: null
    report: null

report =
    name: ''
    order: null
    sample_reports: []
    metric_storage:
        common_for_all_samples_section:
            name: ''
            metrics: []
        sections: []

section =
    name: ''
    metrics: []
    metrics_by_name: {}

metric =
    name: ''
    short_name: ''
    description: ''
    quality: ''
    common: true
    unit: ''


reporting.buildReport = ->
    unless (totalReportData = readJson 'total-report')
        console.log "Error: cannot read #total-report-json"
        return 1

    $('#report_date').html '<p>' + totalReportData.date + '</p>'

    report = totalReportData.report
    metric_storage = report.metric_storage
    metrics_by_name = metric_storage.common_for_all_samples_section.metrics_by_name
    common_records = rec for rec in report.records
    common_records = r for r in common_records when r.metric.name of metrics_by_name
    reporting.buildCommonRecords common_records

    for section of metric_storage.sections
        sample_reports = report.sample_reports
        for rec in sample_reports[0].records when rec.metric.name in (m.name for m in section.metrics)
            console.log rec.metric.name

        columnNames = (rec.metric.name for rec in sample_reports[0].records when rec.metric.name in (m.name for m in section.metrics))
        columnOrder = (recoverOrderFromCookies section.name) or report.order or [0...columnNames.length]

        reporting.buildTotalReport report, section, columnOrder
        plots_html = ""
        for sample_report in sample_reports
            for plot in sample_report.plots
                plots_html += "<img src=\"#{plot}\"/>"
        $('#plot').html plots_html

    return 0

