#!/bin/bash
# Copyright (C) 2011 Ion Torrent Systems, Inc. All Rights Reserved

echo $0

write_html_header ()
{
    local HTML="${RESULTS_DIR}/header"
    if [ -n "$1" ]; then
        HTML="$1"
    fi
    local REFRESHRATE=0
    echo '<html>' > "$HTML"
    print_html_head $REFRESHRATE >> "$HTML"
    echo ' <body>' >> "$HTML"
    print_html_logo >> "$HTML";

    if [ -z "$COV_PAGE_WIDTH" ];then
	echo '  <div id="inner">' >> "$HTML"
    else
	echo "  <div style=\"width:${COV_PAGE_WIDTH}px;margin-left:auto;margin-right:auto;height:100%\">" >> "$HTML"
    fi
    echo "   <h1><center><span style=\"cursor:help\" title=\"Targetted regions='${PLUGIN_TARGETID}'   Target padding=${PLUGIN_PADSIZE}\">Target Coverage Analysis Report</span></center></h1>" >> "$HTML"
}

write_html_footer ()
{
    local HTML="${RESULTS_DIR}/footer"
    if [ -n "$1" ]; then
        HTML="$1"
    fi
    print_html_end_javascript >> "$HTML"
    print_html_footer >> "$HTML"
    echo '  <br/><br/></div>' >> "$HTML"
    echo ' </body></html>' >> "$HTML"
}

display_static_progress ()
{
    local HTML="${RESULTS_DIR}/${HTML_RESULTS}"
    if [ -n "$1" ]; then
        HTML="$1"
    fi
    echo "    <br/><h3 style=\"text-align:center;color:red\">*** Analysis is not complete ***</h3>" >> "$HTML"
    echo "    <a href=\"javascript:document.location.reload();\" ONMOUSEOVER=\"window.status='Refresh'; return true\">" >> "$HTML"
    echo "    <div style=\"text-align:center\">Click here to refresh</div></a>" >> "$HTML"
}

write_json_header ()
{
    local haveBC="false"
    if [ -n "$1" ]; then
	if [ $1 -eq 0 ]; then
	    haveBC="false"
	else
	    haveBC="true"
	fi
    fi
    echo "{" > "$JSON_RESULTS"
    echo "  \"Targetted regions\" : \"$PLUGIN_TARGETID\"," >> "$JSON_RESULTS"
    echo "  \"Target padding\" : \"$PLUGIN_PADSIZE\"," >> "$JSON_RESULTS"
    echo "  \"Examine unique starts\" : \"$PLUGIN_USTARTS\"," >> "$JSON_RESULTS"
    echo -n "  \"barcoded\" : \"$haveBC\"" >> "$JSON_RESULTS"
}

write_json_footer ()
{
    echo -e "\n}" >> "$JSON_RESULTS"
}

write_json_inner ()
{
    local DATADIR=$1
    local DATASET=$2
    local INDENTLEV=$3
    if [ -f "${DATADIR}/${DATASET}/summary.txt" ];then
	echo "," >> "$JSON_RESULTS"
	append_to_json_results "$JSON_RESULTS" "${DATADIR}/${DATASET}/summary.txt" "$DATASET" $INDENTLEV;
    fi
}

