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
        result = null;

    return result;


totalReportData =
    date: null
    reports: []

report =
    name: ''
    order: null
    sample_reports: []


reporting.buildReport = ->
    unless (totalReportData = readJson 'total-report')
        console.log "Error: cannot read #total-report-json"
        return 1

    $('#report_date').html('<p>' + totalReportData.date + '</p>');

    reports = totalReportData.reports
    for report in reports
        report.cornerCell = report.name
    reports[0].cornerCell = 'Sample'

    for report in totalReportData.reports
        sample_reports = report.sample_reports
        columnNames = (record.metric.name for record in sample_reports[0].records)
        columnOrder = (recoverOrderFromCookies report.name) or report.order or [0...columnNames.length]

        reporting.buildTotalReport report, columnOrder

    return 0

