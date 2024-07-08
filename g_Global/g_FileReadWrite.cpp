#include "g_mainwindow.h"
#include <QDebug>
#include <QFile>
#include <QDir>
#include <iostream>
#include <QFileDialog>
#include <QMessageBox>
#include <QStandardPaths>
#include <QSettings>

const char *dir_files="/files";
const char *FileNameAbsent = "/gUntitled";
const char *LISEini="/lisepp.ini";

//QString FFileNameCS;
//QString FFileNameHtml;
QString LISErootPATH;
QString MyDocCompPATH;
QString localPATH;


extern int fontsizeGlobal;
extern int useHighDpiScaling;
FILE *mfopen(const QString& filename, const char* operand);

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void getDPIscaling(void)
{
      MyDocCompPATH = QStandardPaths::standardLocations(QStandardPaths::DocumentsLocation).constFirst();
      QString FN1 = MyDocCompPATH + "/LISEcute";  FN1 += LISEini;
      QSettings myLiseIni1(FN1,QSettings::IniFormat);
      useHighDpiScaling  = myLiseIni1.value("font/scaling",  useHighDpiScaling).toInt();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void getInitialDir(void)
{
    //---------------------------------------------------------  paths begin

//    MyDocCompPATH = QStandardPaths::standardLocations(QStandardPaths::DocumentsLocation).first();
        MyDocCompPATH = QStandardPaths::standardLocations(QStandardPaths::DocumentsLocation).constFirst();
//        qDebug() << MyDocCompPATH;
        QString FileCheck(LISErootPATH); FileCheck += LISEini;
        FILE *fcheck=mfopen(FileCheck,"at");
        int work_in_LISEroot_main = 0;

        if(fcheck) {                    // work in root directory
              fclose(fcheck);
              QSettings myLiseIni0(FileCheck,QSettings::IniFormat);
              work_in_LISEroot_main = myLiseIni0.value("Version/WorkInROOT",0).toInt();
              if(work_in_LISEroot_main) localPATH = LISErootPATH;
              }

        if(work_in_LISEroot_main==0)
              {
              localPATH = MyDocCompPATH;
              localPATH += "/LISEcute";
              }

         //--------------------------------------------------------- lise.ini  begin
          QString FN1=localPATH;  FN1 += LISEini;
          QSettings myLiseIni1(FN1,QSettings::IniFormat);
          fontsizeGlobal     = myLiseIni1.value("font/size",     fontsizeGlobal).toInt();
          useHighDpiScaling  = myLiseIni1.value("font/scaling",  useHighDpiScaling).toInt();
         //--------------------------------------------------------- lise.ini  end

     localPATH += dir_files;
     QDir pathDir(localPATH);
     if(!pathDir.exists()) pathDir.mkdir(localPATH);
    //---------------------------------------------------------  paths end
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::setFileName(QString filename)
{
if(filename.size()>0) FFileName = filename;

QFileInfo info(FFileName);
QString windowName = info.baseName();
if(!windowName.contains(&FileNameAbsent[1])) this->setWindowTitle("Global - " + windowName);
else                                         this->setWindowTitle("Global");


fileNameOut = info.path() + "/" + info.completeBaseName() + ".gloutput";

}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
int MainWindow::checkFileSave(void)
{
int flagVer = QMessageBox::No;
if ( modified ) {
        int flagVer=  QMessageBox::question(this, "Confirmation", "Save Changes?",
                              QMessageBox::Yes | QMessageBox::No | QMessageBox::Cancel,
                              QMessageBox::Yes);

        if(flagVer == QMessageBox::Yes) on_actionSave_triggered();
        }

return flagVer;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_actionOpen_triggered()
{
if(checkFileSave()== QMessageBox::Cancel) return;

QString sfile = QFileDialog::getOpenFileName(this,"Open", FFileName,
            "Global files (*.global);;All files (*.*)");


if(sfile.size()<=0) return;

setFileName(sfile);
readFile(sfile);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::readFile(QString file)
{
    QFile f(file);


FILE *f10 = mfopen(FFileName,"rt");
if(!f10) {
   QMessageBox::information(this, "Open File error","Check file orgigin\n" + FFileName);
   return;
   }
fclose(f10);

QSettings myOpt(FFileName,QSettings::IniFormat);

myOpt.beginGroup("Projectile");
     gAF  = myOpt.value("A",  238).toDouble();
     gZF  = myOpt.value("Z",  92).toDouble();
     gEN0 = myOpt.value("Energy",   430).toDouble();
     gQIN = myOpt.value("Qe",  2).toDouble();
myOpt.endGroup();


myOpt.beginGroup("Target");
   gAT = myOpt.value("A",  63.5).toDouble();
   gZT = myOpt.value("Z",  29).toDouble();
   gDTARGET = myOpt.value("Thick",  100).toDouble();
myOpt.endGroup();


myOpt.beginGroup("Options");
     I_WR     = myOpt.value("Frequency",   2).toInt();
     iOption   = myOpt.value("Option",      0).toInt();
     iCSoutput   = myOpt.value("Output",      2).toInt();
     iLoop   = myOpt.value("Loop",        0).toInt();
     N_Steps  = myOpt.value("N_Steps",    10).toInt();
     Qshow    = myOpt.value("Qshow",       0).toInt();
myOpt.endGroup();

myOpt.beginGroup("Steps");
     DELE   = myOpt.value("E",   30).toDouble();
     DELDT  = myOpt.value("T",  100).toDouble();
     DELZF  = myOpt.value("ZF",   5).toInt();
     DELQ   = myOpt.value("Q",    2).toInt();
     DELZT  = myOpt.value("ZT",   5).toInt();
myOpt.endGroup();
//------------------------------------------------

setPage(true);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_actionSave_triggered()
{
    if(FFileName.size()<=0)
        {
        QMessageBox::information(this, "Save operation error", "File name was not defined\nUse \"Save As\" command");
        //    if ( !strcmp(JustName.c_str(), FileNameAbsent)){SaveAsMenuClick(Sender); return;}
        return;
        }

    FILE *f10=mfopen(FFileName,"wt");
    if(!f10) {
        QMessageBox::information(this, "Open File error","Check file origin\n" + FFileName);
        return;
        }

    fclose(f10);

readPage();

setFileName(FFileName);

QSettings myOpt(FFileName,QSettings::IniFormat);

myOpt.beginGroup("Projectile");
    myOpt.setValue("A",  gAF);
    myOpt.setValue("Z",  gZF);
    myOpt.setValue("Energy",   gEN0);
    myOpt.setValue("Qe",  gQIN);
myOpt.endGroup();


myOpt.beginGroup("Target");
    myOpt.setValue("A",  gAT);
    myOpt.setValue("Z",  gZT);
    myOpt.setValue("Thick",  gDTARGET);
myOpt.endGroup();


myOpt.beginGroup("Options");
    myOpt.setValue("Frequency", I_WR);
    myOpt.setValue("Option",    iOption);
    myOpt.setValue("Output",    iCSoutput);
    myOpt.setValue("Loop",      iLoop);
    myOpt.setValue("N_Steps",   N_Steps);
    myOpt.setValue("Qshow",     Qshow);
myOpt.endGroup();

myOpt.beginGroup("Steps");
    myOpt.setValue("E",    DELE);
    myOpt.setValue("T",    DELDT);
    myOpt.setValue("ZF",   DELZF);
    myOpt.setValue("Q",    DELQ);
    myOpt.setValue("ZT",   DELZT);
myOpt.endGroup();

modified = false;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MainWindow::on_actionSave_As_triggered()
{
    QString sfile = QFileDialog::getSaveFileName(this,  "Save",   FFileName,
                                                 "GLOBAL files (*.global);;All files (*.*)");

    if(sfile.size()<=0) return;

    FFileName = sfile;
    on_actionSave_triggered();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
FILE *mfopen(const QString& filename, const char* operand)
{

#if !defined(__CYGWIN__) && !defined(_WIN32) && !defined(_WIN64)

  FILE *f =   fopen(filename.toStdString().c_str(),operand);

#else

  QString woperand(operand);
  FILE *f = _wfopen(filename.toStdWString().c_str(), woperand.toStdWString().c_str() );

#endif

return f;
}


//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
