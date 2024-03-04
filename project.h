/**************************************************************************/
/* LabWindows/CVI User Interface Resource (UIR) Include File              */
/*                                                                        */
/* WARNING: Do not add to, delete from, or otherwise modify the contents  */
/*          of this include file.                                         */
/**************************************************************************/

#include <userint.h>

#ifdef __cplusplus
    extern "C" {
#endif

     /* Panels and Controls: */

#define  MAIN                             1
#define  MAIN_CHANGE_SWITCH               2       /* control type: binary, callback function: ChangePanel */
#define  MAIN_FREQ_BUTTON                 3       /* control type: command, callback function: OnFrequencyButtonCB */
#define  MAIN_POINTS                      4       /* control type: ring, callback function: (none) */
#define  MAIN_FREQ_PEAK                   5       /* control type: numeric, callback function: (none) */
#define  MAIN_POWER_PEAK                  6       /* control type: numeric, callback function: (none) */
#define  MAIN_GRAPH_RAW_DATA              7       /* control type: graph, callback function: (none) */
#define  MAIN_GRAPH_SPECTRU               8       /* control type: graph, callback function: (none) */
#define  MAIN_GRAPH_RAW_WINDOW            9       /* control type: graph, callback function: (none) */
#define  MAIN_GRAPH_FILTER                10      /* control type: graph, callback function: (none) */
#define  MAIN_GRAPH_SPECTRUM_FILTER       11      /* control type: graph, callback function: (none) */
#define  MAIN_GRAPH_WINDOW                12      /* control type: graph, callback function: (none) */
#define  MAIN_FILTER                      13      /* control type: command, callback function: OnFilterButtonCB */
#define  MAIN_IDC_SECUNDA                 14      /* control type: ring, callback function: (none) */
#define  MAIN_IDC_CUT_FREQ                15      /* control type: numeric, callback function: (none) */
#define  MAIN_IDC_ORDER                   16      /* control type: numeric, callback function: (none) */
#define  MAIN_WINDOW_TYPE                 17      /* control type: ring, callback function: (none) */
#define  MAIN_FILTER_TYPE                 18      /* control type: ring, callback function: (none) */
#define  MAIN_IDC_FPASS                   19      /* control type: numeric, callback function: (none) */

#define  MAIN_PANEL                       2       /* callback function: OnMainPanelCB */
#define  MAIN_PANEL_IDC_GRAPH_FILTER_DATA 2       /* control type: graph, callback function: (none) */
#define  MAIN_PANEL_IDC_GRAPH_RAW_DATA    3       /* control type: graph, callback function: (none) */
#define  MAIN_PANEL_DIM_FEREASTRA         4       /* control type: numeric, callback function: GetVal */
#define  MAIN_PANEL_IDC_TIME_STOP         5       /* control type: numeric, callback function: (none) */
#define  MAIN_PANEL_IDC_TIME_START        6       /* control type: numeric, callback function: (none) */
#define  MAIN_PANEL_FILTER                7       /* control type: ring, callback function: OnFilterClickCB */
#define  MAIN_PANEL_MAX_INDEX_VAL         8       /* control type: numeric, callback function: (none) */
#define  MAIN_PANEL_MIN_INDEX_VAL         9       /* control type: numeric, callback function: (none) */
#define  MAIN_PANEL_ZERO_CROSS_VAL        10      /* control type: numeric, callback function: (none) */
#define  MAIN_PANEL_DISPERSION_VAL        11      /* control type: numeric, callback function: (none) */
#define  MAIN_PANEL_MEDIAN_VAL            12      /* control type: numeric, callback function: (none) */
#define  MAIN_PANEL_MEAN_VAL              13      /* control type: numeric, callback function: (none) */
#define  MAIN_PANEL_MAX_VAL               14      /* control type: numeric, callback function: (none) */
#define  MAIN_PANEL_MIN_VAL               15      /* control type: numeric, callback function: (none) */
#define  MAIN_PANEL_NEXT_BUTTON           16      /* control type: command, callback function: OnNextButtonCB */
#define  MAIN_PANEL_PREV_BUTTON           17      /* control type: command, callback function: OnPrevButtonCB */
#define  MAIN_PANEL_APPLY_BUTTON          18      /* control type: command, callback function: OnApplyButtonCB */
#define  MAIN_PANEL_LOAD_BUTTON           19      /* control type: command, callback function: OnLoadButtonCB */
#define  MAIN_PANEL_COLORNUM2             20      /* control type: numeric, callback function: (none) */
#define  MAIN_PANEL_COLORNUM1             21      /* control type: numeric, callback function: (none) */
#define  MAIN_PANEL_DECORATION_2          22      /* control type: deco, callback function: (none) */
#define  MAIN_PANEL_COMMANDBUTTON         23      /* control type: command, callback function: SaveGraphs */
#define  MAIN_PANEL_BINARYSWITCH2         24      /* control type: binary, callback function: GenDeriv */
#define  MAIN_PANEL_BINARYSWITCH          25      /* control type: binary, callback function: GenAnv */
#define  MAIN_PANEL_GRAPH_HISTOGRAMA      26      /* control type: graph, callback function: (none) */
#define  MAIN_PANEL_CHANGE_SWITCH         27      /* control type: binary, callback function: ChangePanel */
#define  MAIN_PANEL_DECORATION            28      /* control type: deco, callback function: (none) */
#define  MAIN_PANEL_DECORATION_3          29      /* control type: deco, callback function: (none) */


     /* Control Arrays: */

          /* (no control arrays in the resource file) */


     /* Menu Bars, Menus, and Menu Items: */

          /* (no menu bars in the resource file) */


     /* Callback Prototypes: */

int  CVICALLBACK ChangePanel(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK GenAnv(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK GenDeriv(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK GetVal(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK OnApplyButtonCB(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK OnFilterButtonCB(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK OnFilterClickCB(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK OnFrequencyButtonCB(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK OnLoadButtonCB(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK OnMainPanelCB(int panel, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK OnNextButtonCB(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK OnPrevButtonCB(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK SaveGraphs(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);


#ifdef __cplusplus
    }
#endif