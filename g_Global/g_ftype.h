#define Global_version "Version 4.8.1"
#define Global_date    "21-AUG-2026"

// 3.10 Michelle
// 3.11 library


// 4.02 Helmut's changes :  gl_intpoklm.cpp 2017
// 4.03   10/26/20
// Oleg - directory, about, layout, Helmut's change
// 4.04 - version without dll
// 4.05 - version joined with "atima-global"
// 4.6.1  07-01-21  File management
// 4.6.2  07-01-21  Toolbar modification
// 4.6.3  08-01-21  Bugs were fixed
// 4.6.4  08-01-21  Read/Write operations
// 4.6.5  03-02-21  correction in Read/Write operations
// 4.6.6  04-15-21  adaptation of global for Qt 6 --Foster
// 4.6.7  04-22-21  Print and Print Preview actions done --Foster

// 4.6.8  04/27/2021
// revision of paths, creation of arg-file reading utiltiy

// 4.6.9  04/27/2021
// revision of open and save dialogs

// 4.6.10  04/27/2021
// LISErootPATH = QCoreApplication::applicationDirPath();

// 4.6.11  05/21/2021
// corrections in HighDpiScaling attribute at start

// 4.6.12  06/10/2021  adaptation of global for Qt 6 -- Oleg
// 4.6.13  07/12/2021  corrrection for lise.ini path

// 4.6.14  09/07/2021  link on NIM paper
// 4.6.15  12/03/2021  modified for non-latin paths
// 4.6.16  12/03/2021  compatibilty with linux
// 4.6.17  07/13/2022  migrating to Qt63
// 4.6.18  09/07/22  Fixed: bug with mapping in Qt63  (mapped(int) --> mappedInt(int))

// 4.6.19  05/15/23  Global rename of variables to logical names from crazy fortran
// 4.6.20  05/15/23  increase of DELT for thick targets
// 4.6.21  07/01/23  moving to polynom
// 4.6.22  11/06/23  bug correction with sign of a3 in FCORR

// 4.6.23  08/20/26  Modification of project file to avoid McAffee complain on MinGW,
// and use local path for DESTDIR

// 4.7.0   08/20/26  ETACHA-style batch mode from menu.
// Batch file format: Ab Zb Qb Eb At Zt thickness density.
// Results are saved in CSV format.

// 4.7.1   08/20/26  Batch mode writes generated .global and .gloutput files
// to LISEcute/results, forces option #0, and saves summary as *_resGLOBAL.csv.

// 4.7.2   08/20/26  Batch completion message auto-closes only when /q or -q
// command-line option is used.

// 4.7.3   08/21/26  Batch CSV adds iQmax and dQ(calc-set), where
// dQ(calc-set) = Qcalc - Qb.

// 4.7.4   08/21/26  Batch CSV renames Result to Status, reports iQmax as
// the charge state with maximum fraction, initializes result buffers, and
// marks rows without physical Qcalc output as Status 0.

// 4.7.5   08/21/26  Batch mode treats ETACHA Qb as charge state and converts
// it to GLOBAL's Z-q input before each calculation.

// 4.7.6   08/21/26  Batch CSV fraction output is limited to Fq0 through Fq9.

// 4.8.0   08/21/26  Numerical refinements proposed by Masahiro Yoshimoto,
// RIKEN BigRIPS team <masahiro.yoshimoto@riken.jp>:
// EDIFF divisor changed from 200 to 2000 and integration coefficients changed
// from {1,10,100,1000,10000,100000} to {1,1,10,100,1000,10000}.

// 4.8.1   08/21/26  ETACHA-compatible command-line batch start:
// Global.exe batch-file -r or /r loads the batch file and immediately runs it.
