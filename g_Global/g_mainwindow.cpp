#include "g_mainwindow.h"
#include "ui_g_mainwindow.h"

#include <QFile>
#include <QTextStream>
#include <QMessageBox>
#include <QPrinter>
#include <QPrintDialog>
#include <QDesktopServices>
#include <QDebug>
#include <QSignalMapper>
#include <QToolButton>
#include <QPainter>
#include <QPrintPreviewDialog>
#include <QDesktopServices>
#include <QUrl>
#include <QDir>
#include <QFileDialog>
#include <QFileInfo>
#include <QLocale>
#include <QRegularExpression>
#include <QVector>
#include <QTimer>


#include "L_Init/Constant.h"
#include "g_ftype.h"
#include "g_about.h"
#include "L_Atima/global/globallib.h"

#define aem_mg      1.6605402e-21
#define U_FC        aem_mg


extern void ELEMENT(double &Z, double &A, char *CZ, int IOPT, int &IRC );
extern int RunGlobalLocal(double *Obmen, char *filename, bool option_read_data_file,
                        int fast,  char *Globalversion,
                        double gAF, double gZF, double gAT, double gZT,
                        double gDTARGET, double gEN0, int  gQIN,
                        int   I_WR, int   iOption, int   iCSoutput,
                        int   iLoop, int N_Steps, int  Qshow,
                        int    DELZF, double DELE,
                        int    DELQ, int DELZT, double DELDT);

int runGlobal(char *filename, bool option_read_data_file, const char *Globalversion,
               double gAF, double gZF, double gAT, double gZT,
               double gDTARGET, double gEN0, int gQIN, int I_WR, int iOption,
              int iCSoutput, int iLoop, int N_Steps, int Qshow,
              int DELZF, double DELE, int DELQ, int DELZT, double DELDT,
              double *ObmenOut = nullptr);

double D1;

extern QString LISErootPATH;
extern QString localPATH;
extern const char *FileNameAbsent;
extern QString FileArg;
extern bool AutoCloseBatchMessage;
extern bool AutoRunBatchFile;

class batchReactionSet
{
public:
    bool initFromString(const QString& line);

    int Ab = 0;
    int Zb = 0;
    int Qb = 0;
    double Eb = 0.0;
    int At = 0;
    int Zt = 0;
    double thick = 0.0;
    double density = 0.0;
};

bool batchReactionSet::initFromString(const QString& line)
{
#if QT_VERSION >= QT_VERSION_CHECK(5, 14, 0)
    QStringList parts = line.split(QRegularExpression("[,\\s]+"), Qt::SkipEmptyParts);
#else
    QStringList parts = line.split(QRegularExpression("[,\\s]+"), QString::SkipEmptyParts);
#endif

    if(parts.size() < 8) return false;

    bool ok = true;
    int i = 0;

    Ab      = parts[i++].toInt(&ok);     if(!ok || Ab      <= 0) return false;
    Zb      = parts[i++].toInt(&ok);     if(!ok || Zb      <= 0) return false;
    Qb      = parts[i++].toInt(&ok);     if(!ok || Qb      < 0 || Qb > Zb) return false;
    Eb      = parts[i++].toDouble(&ok);  if(!ok || Eb      <= 0.0) return false;
    At      = parts[i++].toInt(&ok);     if(!ok || At      <= 0) return false;
    Zt      = parts[i++].toInt(&ok);     if(!ok || Zt      <= 0) return false;
    thick   = parts[i++].toDouble(&ok);  if(!ok || thick   <= 0.0) return false;
    density = parts[i++].toDouble(&ok);  if(!ok || density <= 0.0) return false;

    return true;
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ui->statusBar->hide();

    gAF = 238.;
    gZF = 92.;
    gAT = 63.5;
    gZT = 29.;
    gDTARGET = 100.;
    gEN0 = 430.;
    gQIN = 2;


    I_WR  =  0;        // frequency
    iOption = opt_NormIni;
    iCSoutput = 2;
    iLoop = 0;
    N_Steps = 10;
    Qshow = 0;

    DELZF = 5;
    DELE = 30.;
    DELQ =2;
    DELZT = 5;
    DELDT = 100;

//----------------------------------------------------------------------
    QSignalMapper *signalMapper1 = new QSignalMapper(this);


    connect(ui->radioButton_1,SIGNAL(clicked()),signalMapper1, SLOT(map()));    signalMapper1->setMapping(ui->radioButton_1,0);
    connect(ui->radioButton_2,SIGNAL(clicked()),signalMapper1, SLOT(map()));    signalMapper1->setMapping(ui->radioButton_2,1);
    connect(ui->radioButton_3,SIGNAL(clicked()),signalMapper1, SLOT(map()));    signalMapper1->setMapping(ui->radioButton_3,2);
    connect(ui->radioButton_4,SIGNAL(clicked()),signalMapper1, SLOT(map()));    signalMapper1->setMapping(ui->radioButton_4,3);
    connect(ui->radioButton_5,SIGNAL(clicked()),signalMapper1, SLOT(map()));    signalMapper1->setMapping(ui->radioButton_5,4);
    connect(ui->radioButton_6,SIGNAL(clicked()),signalMapper1, SLOT(map()));    signalMapper1->setMapping(ui->radioButton_6,5);

    connect(signalMapper1, SIGNAL(mappedInt(int)),this, SLOT(loopChanged(int)));
//----------------------------------------------------------------------
    QSignalMapper *signalMapper2 = new QSignalMapper(this);
    connect(ui->radioButton_7,SIGNAL(clicked()),signalMapper2, SLOT(map()));    signalMapper2->setMapping(ui->radioButton_7,0);
    connect(ui->radioButton_8,SIGNAL(clicked()),signalMapper2, SLOT(map()));    signalMapper2->setMapping(ui->radioButton_8,1);
    connect(ui->radioButton_9,SIGNAL(clicked()),signalMapper2, SLOT(map()));    signalMapper2->setMapping(ui->radioButton_9,2);
    connect(ui->radioButton_10,SIGNAL(clicked()),signalMapper2, SLOT(map()));   signalMapper2->setMapping(ui->radioButton_10,3);
    connect(ui->radioButton_11,SIGNAL(clicked()),signalMapper2, SLOT(map()));   signalMapper2->setMapping(ui->radioButton_11,4);

    connect(signalMapper2,SIGNAL(mappedInt(int)),this, SLOT(optionsChanged(int)));
//----------------------------------------------------------------------

    ZPedit_permit=true;
    ZTedit_permit=true;

    modified=false;

    FFileName  = FileArg.size()>0? FileArg : localPATH + FileNameAbsent;
    setFileName(FFileName);
    if(FileArg.size()>0 && !AutoRunBatchFile) readFile(FFileName);

    setPage();
    if(FileArg.size()>0 && AutoRunBatchFile)
        QTimer::singleShot(0, this, [this]() { runBatchFile(FFileName); });
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
MainWindow::~MainWindow(){    delete ui;}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::setPage(bool keepA)
{
double   gAFs=gAF, gATs=gAT;

    ZPedit_permit = false;
    ui->proj_A->setText(QString::number(gAF,'f',0));
    ui->proj_Z->setText(QString::number(gZF,'f',0));
    ui->proj_ZQ->setText(QString::number(gQIN));
    ZPedit_permit = true;

    ZTedit_permit = false;
    ui->targ_A->setText(QString::number(gAT,'f',2));//,'f',0));
    ui->targ_Z->setText(QString::number(gZT,'f',0));//,'f',0));
    ui->initEnergy->setText(QString::number(gEN0));//,'f',0));
    ui->thickness->setText(QString::number(gDTARGET));//,'f',0));
    ZTedit_permit = true;

    on_proj_Z_textEdited(QString::number(gZF));
    on_targ_Z_textEdited(QString::number(gZT));

if(keepA)
    {
    gAF=gAFs; gAT=gATs;
    ui->targ_A->setText(QString::number(gAT,'f',2));//,'f',0));
    ui->proj_A->setText(QString::number(gAF,'f',0));
    }

    ui->num_steps->setText(QString::number(N_Steps));
    ui->plotQ->setText(QString::number(Qshow));
    switch(iLoop) {
    case 0:
        ui->radioButton_1->setChecked(true);
        break;
    case 1:
        ui->radioButton_2->setChecked(true);
        break;
    case 2:
        ui->radioButton_3->setChecked(true);
        break;
    case 3:
        ui->radioButton_4->setChecked(true);
        break;
    case 4:
        ui->radioButton_5->setChecked(true);
        break;
    case 5:
        ui->radioButton_6->setChecked(true);
        break;
    }
    switch(iOption) {
    case 0:
        ui->radioButton_7->setChecked(true);
        break;
    case 1:
        ui->radioButton_8->setChecked(true);
        break;
    case 2:
        ui->radioButton_9->setChecked(true);
        break;
    case 3:
        ui->radioButton_10->setChecked(true);
        break;
    case 4:
        ui->radioButton_11->setChecked(true);
        break;
    }

switch(I_WR) {
    case 0:        ui->radioButton_12->setChecked(true);        break;
    case 1:        ui->radioButton_13->setChecked(true);        break;
    case 2:        ui->radioButton_14->setChecked(true);        break;
    case 3:        ui->radioButton_15->setChecked(true);        break;
    }

    ui->cb_CSoutput->setChecked(iCSoutput==3 ? true : false);
    ui->dZ_proj->setText(QString::number(DELZF));
    ui->energy_step->setText(QString::number(DELE));
    ui->dQstep->setText(QString::number(DELQ));
    ui->dZ_target->setText(QString::number(DELZT));
    ui->thick_step->setText(QString::number(DELDT));
    loopChanged(iLoop);
    optionsChanged(iOption);
    MakeTarget();
    modified = false;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::readPage()
{
QString cc[8] = {"Warning: wrong initial settings",
                 "Energy should be more than 30 MeV/u\n\r"
                 "and less than 2000 MeV/u!",

                 "ZP < 29! GLOBAL is developed for ZP > 28!\n"
                 "Check result",

                 "Check Q-state!",
                 "One of delta-steps is wrong",
                 "Thickness should be positive!",
                 "Number of Steps > 0 and < 100!",
                 "Q-array starts from 0 up to 25!"
                 };


frequencyChanged();

iCSoutput = (ui->cb_CSoutput->isChecked() == true ? 3 : 2);
gEN0 = ui->initEnergy->text().toDouble();
if(gEN0 < 30 || gEN0 > 2000)
    {
    QMessageBox MB1;
    MB1.setWindowTitle(cc[0]);
    MB1.setText(cc[1]);
    MB1.setIcon(QMessageBox::Warning);
    MB1.exec();
    return;
}
if(gZF < 29) {
    QMessageBox MB1;
    MB1.setWindowTitle(cc[0]);
    MB1.setText(cc[2]);
    MB1.setIcon(QMessageBox::Warning);
    MB1.exec();
    return;
}
gQIN = ui->proj_ZQ->text().toInt();
if(gQIN<0 || gQIN>gZF || gQIN > 27){
    QMessageBox MB1;
    MB1.setWindowTitle(cc[0]);
    MB1.setText(cc[3]);
    MB1.setIcon(QMessageBox::Warning);
    MB1.exec();
    return;
}
gAF = ui->proj_A->text().toDouble();
DELZF = ui->dZ_proj->text().toInt();
DELE = ui->energy_step->text().toDouble();
DELQ = ui->dQstep->text().toInt();
DELZT = ui->dZ_target->text().toInt();
DELDT = ui->thick_step->text().toDouble();
N_Steps = ui->num_steps->text().toInt();

bool flagStep = false;
if(DELZF <= 0) {DELZF=5;  flagStep=true;};
if(DELE  <= 0) {DELE=30;  flagStep=true;};
if(DELQ  <= 0) {DELQ=5;   flagStep=true;};
if(DELZT <= 0) {DELZT=5;  flagStep=true;};
if(DELDT <= 0) {DELDT=90; flagStep=true;};

if(flagStep && iLoop!=0){
    QMessageBox MB1;
    MB1.setWindowTitle(cc[0]);
    MB1.setText(cc[4]);
    MB1.setIcon(QMessageBox::Warning);
    MB1.exec();
}
if(gDTARGET<=0){
    QMessageBox MB1;
    MB1.setWindowTitle(cc[0]);
    MB1.setText(cc[5]);
    MB1.setIcon(QMessageBox::Warning);
    MB1.exec();
}
if( (N_Steps <= 0 || N_Steps > 100) &&
    (iLoop==2  || iLoop==5) )
{
    N_Steps = 5;
    QMessageBox MB1;
    MB1.setWindowTitle(cc[0]);
    MB1.setText(cc[6]);
    MB1.setIcon(QMessageBox::Warning);
    MB1.exec();
}
Qshow = ui->plotQ->text().toInt();
if( (Qshow < 0 || Qshow > 25))    {
        Qshow=0;
        QMessageBox MB1;
        MB1.setWindowTitle(cc[0]);
        MB1.setText(cc[7]);
        MB1.setIcon(QMessageBox::Warning);
        MB1.exec();
    }
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::MakeTarget()
{
    double TTHICK_A = gDTARGET /(gAT*U_FC);     // atoms/cm2;
    ui->atoms->setText(QString::number(TTHICK_A,'e',4));
    modified = true; ui->textEdit->setText("");

}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_actionExecute_triggered()
{
if(fileNameOut.size()<=0)    {
    QMessageBox::information(this,"Error!", "output file name is empty");
    return;
    }

    char tGlobal_version[20];
    strcpy(tGlobal_version, Global_version);
    QByteArray FNO = fileNameOut.toLocal8Bit();
    char *fno = FNO.data();
    ui->label_ready->setText("<h3 style=\"color: #f00;\">One moment please...</h3>");
    QCoreApplication::processEvents();
    readPage();

result = runGlobal (fno, false, tGlobal_version,
                   gAF,  gZF,  gAT,  gZT,
                   gDTARGET,  gEN0,   gQIN,
                     I_WR,    iOption,    iCSoutput,
                     iLoop,  N_Steps,   Qshow,
                      DELZF,  DELE,
                      DELQ,  DELZT,  DELDT);

QString content;

if(result==-697)
    {
    QMessageBox::information(this,"File open error", FFileName);
    }
else {
    ui->label_ready->setText("<h3 style=\"color: #00f;\">Ready</h3>");
    QFile resFile(fileNameOut);
    resFile.open(QIODevice::ReadOnly);
    QTextStream stream(&resFile);
    content = stream.readAll();
    resFile.close();
    }

ui->textEdit->setText(content);
ui->textEdit->setReadOnly(true);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_actionBatch_mode_triggered()
{
    QString openFrom = BatchFileName.isEmpty() ? localPATH : BatchFileName;
    QString fileName = QFileDialog::getOpenFileName(
                this,
                tr("Open GLOBAL batch file"),
                openFrom,
                tr("GLOBAL batch files (*.bglobal *.betacha);;All files (*.*)"));

    if(fileName.isEmpty()) return;

    BatchFileName = fileName;
    runBatchFile(fileName);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
bool MainWindow::runBatchFile(const QString &fileName)
{
    QFile file(fileName);
    if(!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        QMessageBox::warning(this, tr("File error"),
                             tr("Cannot open batch file:\n%1").arg(fileName));
        return false;
    }

    QVector<batchReactionSet> reactions;
    QTextStream in(&file);
    int lineCount = 0;
    int commentCount = 0;
    while(!in.atEnd()) {
        QString line = in.readLine().trimmed();
        lineCount++;

        if(line.isEmpty() || line.startsWith('!') || line.startsWith(';')) {
            commentCount++;
            continue;
        }

        batchReactionSet r;
        if(r.initFromString(line)) reactions.append(r);
    }
    file.close();

    if(reactions.isEmpty()) {
        QMessageBox::information(this, tr("No data"),
                                 tr("No valid data lines were found in:\n%1").arg(fileName));
        return false;
    }

    QFileInfo batchInfo(fileName);
    QDir batchDir = batchInfo.absoluteDir();
    QString resultPath = batchDir.filePath(batchInfo.completeBaseName() + "_resGLOBAL.csv");
    QDir resultsDir(QFileInfo(localPATH).absoluteDir().filePath("results"));
    if(!resultsDir.exists()) resultsDir.mkpath(".");

    QFile outFile(resultPath);
    if(!outFile.open(QIODevice::WriteOnly | QIODevice::Text)) {
        QMessageBox::warning(this, tr("File error"),
                             tr("Cannot open results file for writing:\n%1").arg(resultPath));
        return false;
    }

    QTextStream ts(&outFile);
    ts.setLocale(QLocale::C);
    ts.setRealNumberNotation(QTextStream::FixedNotation);
    ts.setRealNumberPrecision(6);
    auto csvText = [](QString value) {
        value.replace('"', "\"\"");
        return "\"" + value + "\"";
    };

    ts << "Ab,Zb,Qb,Eb,At,Zt,thick,density,"
       << "Status,Qcalc,dQ,iQmax,dQ(calc-set),EqThick,EnergyOut,";
    for(int q = 0; q <= 9; ++q) ts << "Fq" << q << ',';
    ts << "globalFile,outputFile\n";

    QString oldFFileName = FFileName;
    QString oldFileNameOut = fileNameOut;
    bool oldModified = modified;
    double oldAF = gAF;
    double oldZF = gZF;
    double oldAT = gAT;
    double oldZT = gZT;
    double oldDTARGET = gDTARGET;
    double oldEN0 = gEN0;
    int oldQIN = gQIN;
    int oldOption = iOption;
    int oldLoop = iLoop;

    ui->label_ready->setText("<h3 style=\"color: #f00;\">Batch mode...</h3>");
    ui->menuBar->setEnabled(false);
    QCoreApplication::processEvents();

    int completedSets = 0;
    for(int i = 0; i < reactions.size(); ++i) {
        const batchReactionSet &r = reactions[i];

        gAF = r.Ab;
        gZF = r.Zb;
        gQIN = r.Zb - r.Qb;
        gEN0 = r.Eb;
        gAT = r.At;
        gZT = r.Zt;
        gDTARGET = r.thick;
        iOption = opt_NormIni;
        iLoop = 0;

        setPage(true);

        QString beamEl = ui->proj_element->text();
        beamEl.remove(' ');
        QString targEl = ui->targ_element->text();
        targEl.remove(' ');

        QString rowBase = QString::asprintf(
                    "%d%s%d_e%.1f_%d%s_t%.1f",
                    r.Ab, beamEl.toUtf8().constData(), r.Qb, r.Eb,
                    r.At, targEl.toUtf8().constData(), r.thick);

        QString globalFile = resultsDir.filePath(rowBase + ".global");
        QString outputFile = resultsDir.filePath(rowBase + ".gloutput");

        if(!saveGlobalFile(globalFile)) {
            QMessageBox::warning(this, tr("File error"),
                                 tr("Cannot write generated GLOBAL file:\n%1").arg(globalFile));
            continue;
        }

        char tGlobal_version[20];
        strcpy(tGlobal_version, Global_version);
        QByteArray FNO = outputFile.toLocal8Bit();
        double Obmen[IMAX+4];
        for(int j = 0; j < IMAX+4; ++j) Obmen[j] = 0.0;

        int batchResult = runGlobal(FNO.data(), false, tGlobal_version,
                                    gAF, gZF, gAT, gZT,
                                    gDTARGET, gEN0, gQIN,
                                    I_WR, opt_NormIni, iCSoutput,
                                    0, N_Steps, Qshow,
                                    DELZF, DELE,
                                    DELQ, DELZT, DELDT,
                                    Obmen);

        double qCalc = Obmen[IMAX];
        int csvStatus = batchResult;
        if(csvStatus == 1 && qCalc <= 0.0) csvStatus = 0;

        int iQmax = 0;
        if(qCalc > 0.0) {
            int electronIndexMax = 0;
            for(int q = 1; q < IMAX; ++q)
                if(Obmen[q] > Obmen[electronIndexMax]) electronIndexMax = q;
            iQmax = r.Zb - electronIndexMax;
        }

        double dQcalcSet = qCalc - r.Qb;

        ts << r.Ab << ',' << r.Zb << ',' << r.Qb << ',' << r.Eb << ','
           << r.At << ',' << r.Zt << ',' << r.thick << ',' << r.density << ','
           << csvStatus << ',' << qCalc << ',' << Obmen[IMAX+1] << ','
           << iQmax << ',' << dQcalcSet << ','
           << Obmen[IMAX+2] << ',' << Obmen[IMAX+3] << ',';
        for(int q = 0; q <= 9; ++q) ts << Obmen[q] << ',';
        ts << csvText(globalFile) << ',' << csvText(outputFile) << '\n';

        completedSets++;
        ui->label_ready->setText(tr("<h3 style=\"color: #f00;\">Batch %1/%2...</h3>")
                                 .arg(i + 1).arg(reactions.size()));
        QCoreApplication::processEvents();
    }

    outFile.close();
    ui->menuBar->setEnabled(true);
    FFileName = oldFFileName;
    fileNameOut = oldFileNameOut;
    gAF = oldAF;
    gZF = oldZF;
    gAT = oldAT;
    gZT = oldZT;
    gDTARGET = oldDTARGET;
    gEN0 = oldEN0;
    gQIN = oldQIN;
    iOption = oldOption;
    iLoop = oldLoop;
    modified = oldModified;
    setPage(true);
    ui->label_ready->setText("<h3 style=\"color: #00f;\">Ready</h3>");

    QMessageBox *messageBox = new QMessageBox(
                QMessageBox::Information,
                tr("Batch completed"),
                tr("Completed %1 calculations from %2 lines (%3 comments/empty).\nResults saved to:\n%4")
                    .arg(completedSets)
                    .arg(lineCount)
                    .arg(commentCount)
                    .arg(QDir::toNativeSeparators(resultPath)),
                QMessageBox::Ok,
                this);
    messageBox->setAttribute(Qt::WA_DeleteOnClose);

    if(AutoCloseBatchMessage) {
        QTimer::singleShot(5000, messageBox, &QMessageBox::accept);
        messageBox->open();
    } else {
        messageBox->exec();
    }

    QFile resFile(resultPath);
    if(resFile.open(QIODevice::ReadOnly | QIODevice::Text)) {
        QTextStream stream(&resFile);
        ui->textEdit->setText(stream.readAll());
        ui->textEdit->setReadOnly(true);
        resFile.close();
    }

    return true;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_proj_Z_textEdited(const QString &arg1)
{
    if(!ZPedit_permit) return;
    if(arg1.length() == 0) return;
    gZF = arg1.toInt();
    strcpy(CZ,"**");
    int IRC1;
    if(gZF>0 && gZF <=97)
        {
        ELEMENT(gZF, gAF, CZ, 3, IRC1);
        if(IRC1 == -1) {
            QMessageBox MB1;
            MB1.setWindowTitle("Error!");
            MB1.setText("Z <-> Element conversion.");
            MB1.setIcon(QMessageBox::Warning);
            MB1.exec();
            }
        }
    ZPedit_permit = false;

    ui->proj_element->setText(QString::fromLocal8Bit(CZ));
    ui->proj_A->setText(QString::number(gAF,'f',0));
    ZPedit_permit = true;
    CheckQState();
    modified = true; ui->textEdit->setText("");

}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_targ_Z_textEdited(const QString &arg1)
{
    if(!ZTedit_permit) return;
    if(arg1.length() == 0) return;
    gZT = arg1.toInt();
    strcpy(CZ,"**");
    int IRC1;
    if(gZT>0 && gZT <=97) {
        ELEMENT(gZT, gAT, CZ, 1, IRC1);
        if(IRC1 == -1) {
            QMessageBox MB1;
            MB1.setWindowTitle("Error!");
            MB1.setText("Z <-> Element conversion.");
            MB1.setIcon(QMessageBox::Warning);
            MB1.exec();
        }
    }
    ui->targ_element->setText(QString::fromLocal8Bit(CZ));
    ui->targ_A->setText(QString::number(gAT,'f',2));
    ZTedit_permit = true;
    MakeTarget();
    modified = true; ui->textEdit->setText("");
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_proj_element_textEdited(const QString &arg1)
{
   if(!ZPedit_permit) return;

   if(arg1.length() == 0) return;
   strcpy(CZ,arg1.toLocal8Bit());
   if(strlen(CZ)==1) strcat(CZ," ");
   int IRC1;

   ELEMENT(gZF, gAF, CZ, 4, IRC1 );
   ZPedit_permit=false;
   if(IRC1==0){
       ui->proj_Z->setText(QString::number(gZF,'f',0));
       ui->proj_A->setText(QString::number(gAF,'f',0));
   } else {
       ui->proj_Z->setText("**");
       ui->proj_A->setText("**");
   }

   ZPedit_permit = true;
   CheckQState();
   modified = true; ui->textEdit->setText("");

}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_action_Exit_triggered()
{
if(checkFileSave()== QMessageBox::Cancel) return;
exit(1);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_actionPrint_triggered()
{
    QPrinter printer(QPrinter::HighResolution);
    QPrintDialog *printDialog = new QPrintDialog(&printer, this);
    printDialog->setWindowTitle(tr("Print Widget"));

    if (printDialog->exec() != QDialog::Accepted)
        return;
    QPainter painter;
    painter.begin(&printer);



    double xscale = printer.pageRect(QPrinter::DevicePixel).width()/double(this->width());
    double yscale = printer.pageRect(QPrinter::DevicePixel).height()/double(this->height());
    double scale = 0.9 * qMin(xscale, yscale);
    painter.translate(printer.paperRect(QPrinter::DevicePixel).x() + printer.pageRect(QPrinter::DevicePixel).width()/2,
                      printer.paperRect(QPrinter::DevicePixel).y() + printer.pageRect(QPrinter::DevicePixel).height()/2);


    painter.scale(scale, scale);
    painter.translate(-width()/2, -height()/2);

    this->render(&painter);
    painter.end();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_actionWeb_Documentation_triggered()
{
    QDesktopServices::openUrl(QUrl("http://lise.nscl.msu.edu/6_3/lise++_6_3.pdf#page=2"));
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_thickness_textEdited(const QString &arg1)
{
    gDTARGET = (arg1.length() > 0 ? arg1.toDouble() : 0);
    modified = true;
    MakeTarget();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_initEnergy_textEdited(const QString &arg1)
{
    gEN0 = arg1.toDouble();
    modified = true; ui->textEdit->setText("");

}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::CheckQState()
{
gQIN = ui->proj_ZQ->text().toInt();
if(gQIN<0 || gQIN>gZF)
        {
        gQIN = 0;
        ui->proj_ZQ->setText(QString::number(gQIN));
        }

}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_targ_element_textEdited(const QString &arg1)
{
    if(!ZTedit_permit) return;

    if(arg1.length() == 0) return;
    strcpy(CZ,arg1.toLocal8Bit());
    if(strlen(CZ)==1) strcat(CZ," ");
    int IRC1;

    ELEMENT(gZT, gAT, CZ, 2, IRC1 );
    ZTedit_permit=false;
    if(IRC1==0){
        ui->targ_Z->setText(QString::number(gZT,'f',0));
        ui->targ_A->setText(QString::number(gAT,'f',2));
    } else {
        ui->targ_Z->setText("**");
        ui->targ_A->setText("**");
    }

    ZTedit_permit = true;
    modified = true; ui->textEdit->setText("");
    MakeTarget();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_targ_A_textEdited(const QString &arg1)
{
    if(!ZTedit_permit) return;
    gAT = arg1.toDouble();
    MakeTarget();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
 void MainWindow::loopChanged(int i)
 {
   //  qDebug() << "loopChanged";
   //  qDebug() << i;
     iLoop = i;
     ui->dZ_proj->setVisible(iLoop==1);
     ui->label_dZ->setVisible(iLoop==1);

     ui->label_Estep->setVisible(iLoop ==2);
     ui->energy_step->setVisible(iLoop ==2);

     ui->label_dQ->setVisible(iLoop ==3);
     ui->dQstep->setVisible(iLoop ==3);

     ui->label_dZtarget->setVisible(iLoop==4);
     ui->dZ_target->setVisible(iLoop==4);

     ui->label_thickStep->setVisible(iLoop==5);
     ui->thick_step->setVisible(iLoop==5);

     ui->label_nSteps->setVisible(iLoop==5 || iLoop==2);
     ui->num_steps->setVisible(iLoop==5 || iLoop==2);

     if(iLoop==2)
         {
         bool FlagEfinal = (iOption==opt_EqEnd || iOption==opt_NormEnd);
         if(FlagEfinal)
             {
             iOption--;
                   if(iOption==0){ui->radioButton_7->setChecked(true);  optionsChanged(0);}
             else  if(iOption==1){ui->radioButton_8->setChecked(true);  optionsChanged(1);}
             else  if(iOption==2){ui->radioButton_9->setChecked(true);  optionsChanged(2);}
             else  if(iOption==3){ui->radioButton_10->setChecked(true); optionsChanged(3);}
             else  if(iOption==4){ui->radioButton_11->setChecked(true); optionsChanged(4);}
             }
        }

 modified = true; ui->textEdit->setText("");
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
 void MainWindow::optionsChanged(int i)
 {
     iOption = i;
     ui->freq_groupBox->setEnabled(iOption==opt_Evolution);

     bool thick_flag = (iOption!=opt_EqIni  && iOption!=opt_EqEnd);
     bool FlagEfinal = (iOption==opt_EqEnd || iOption==opt_NormEnd);

     ui->thickness->setEnabled(thick_flag);
     ui->thick_step->setEnabled(thick_flag);
     QString fe = (FlagEfinal ? "Final" : "Initial");
     ui->label_energy->setText("<b>" + fe+"</b> &nbsp;  Energy");

//------
     if(iLoop==2)
         {
         bool FlagEfinal = (iOption==opt_EqEnd || iOption==opt_NormEnd);
         if(FlagEfinal)
             {
             iOption--;
                   if(iOption==0){ui->radioButton_7->setChecked(true);  optionsChanged(0);}
             else  if(iOption==1){ui->radioButton_8->setChecked(true);  optionsChanged(1);}
             else  if(iOption==2){ui->radioButton_9->setChecked(true);  optionsChanged(2);}
             else  if(iOption==3){ui->radioButton_10->setChecked(true); optionsChanged(3);}
             else  if(iOption==4){ui->radioButton_11->setChecked(true); optionsChanged(4);}
             }
        }
//------

     modified = true; ui->textEdit->setText("");

 }
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::frequencyChanged()
{
     if(ui->radioButton_12->isChecked()) I_WR=0;
else if(ui->radioButton_13->isChecked()) I_WR=1;
else if(ui->radioButton_14->isChecked()) I_WR=2;
else                                     I_WR=3;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
 void MainWindow::on_actionAbout_triggered()
{
    About *about_page = new About;
    about_page->setWindowFlags(Qt::CustomizeWindowHint |
                               Qt::WindowTitleHint | Qt::WindowCloseButtonHint);
    about_page->show();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
 void MainWindow::keyPressEvent(QKeyEvent *e)
{
if(e->key()==Qt::Key_Enter || e->key()==Qt::Key_Return){
    on_actionExecute_triggered();
    }
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
 int runGlobal(char *filename, bool option_read_data_file, const char *Globalversion,
                double gAF, double gZF, double gAT, double gZT,
                double gDTARGET, double gEN0, int gQIN, int I_WR, int iOption,
               int iCSoutput, int iLoop, int N_Steps, int Qshow,
               int DELZF, double DELE, int DELQ, int DELZT, double DELDT,
               double *ObmenOut)
 {
     double Obmen[IMAX+4];    //Qmean, dQ, EquilibThick, Eout
     for(int i=0; i<IMAX+4; i++) Obmen[i] = 0.;
     int result = RunGlobalLocal(Obmen,filename, option_read_data_file,
                             0,(char*)Globalversion,
                             gAF, gZF, gAT, gZT,
                             gDTARGET, gEN0, gQIN,
                             I_WR, iOption, iCSoutput,
                             iLoop, N_Steps, Qshow,
                             DELZF, DELE,
                             DELQ,  DELZT, DELDT);
     if(ObmenOut)
       for(int i=0; i<IMAX+4; i++) ObmenOut[i] = Obmen[i];
     return result;
 }
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_actionPrint_Preview_triggered()
{
#ifndef QT_NO_PRINTER
    QPrinter printer(QPrinter::HighResolution);
    QPrintPreviewDialog preview(&printer, this);
    int width = 1.2 * printer.pageRect(QPrinter::Point).width();
    int height = 1.2 * printer.pageRect(QPrinter::Point).height();
    preview.setMinimumSize(width,height);
    preview.setWindowFlags ( Qt::Window );
    setWindowTitle(tr("Print Preview"));
    connect(&preview, SIGNAL(paintRequested(QPrinter*)),this, SLOT(printPreview(QPrinter*)));
    QPainter painter;
    preview.exec();
#endif
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::printPreview(QPrinter* printer)
{
    QPainter painter(printer);
    painter.setRenderHints(QPainter::Antialiasing |
                           QPainter::TextAntialiasing |
                           QPainter::SmoothPixmapTransform, true);
    double xscale = printer->pageRect(QPrinter::DevicePixel).width()/double(this->width());
    double yscale = printer->pageRect(QPrinter::DevicePixel).height()/double(this->height());
    double scale = 0.9 * qMin(xscale, yscale);
    painter.translate(printer->paperRect(QPrinter::DevicePixel).x() + printer->pageRect(QPrinter::DevicePixel).width()/2,
                      printer->paperRect(QPrinter::DevicePixel).y() + printer->pageRect(QPrinter::DevicePixel).height()/2);


    painter.scale(scale, scale);
    painter.translate(-width()/2, -height()/2);

    this->render(&painter);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_actionNIM_triggered()
{
 {    QDesktopServices::openUrl(QUrl("http://lise.nscl.msu.edu/doc/charge-global.pdf"));}
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
