#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QKeyEvent>
#include <QPrinter>
//#include "g_global.h"
//------------------------------------------------
namespace Ui {
class MainWindow;
}
//------------------------------------------------
class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private:

    bool ZPedit_permit;
    bool ZTedit_permit;

    int result;
    bool modified;

    QString FFileName;
    QString fileNameOut;

    char CZ[20];
    double gAF;
    double gZF;
    double gAT;
    double gZT;
    double gDTARGET;
    double gEN0;
    int    gQIN;

    int   I_WR;
    int   iOption;// = opt_NormIni;
    int   iCSoutput;
    int   iLoop;
    int   N_Steps;
    int   Qshow;

    double DELE;
    double DELDT;
    int    DELZF;
    int    DELQ;
    int    DELZT;


    void readPage();
    void setPage(bool keepA=false);
    void CheckQState();
//    void ELEMENT(double &Z, double &A, char *CZ, int IOPT, int &IRC );
    void MakeTarget();
    void setFileName(QString filename);
//    void getInitialDir();
    int checkFileSave();
    void frequencyChanged();
    void readFile(QString file);


protected:
    void keyPressEvent(QKeyEvent *e);

//------------------------------------------------

//public slots:
private slots:
    void loopChanged(int);
    void optionsChanged(int);


private slots:
    void on_proj_Z_textEdited(const QString &arg1);
    void on_targ_Z_textEdited(const QString &arg1);
    void on_proj_element_textEdited(const QString &arg1);

    void on_thickness_textEdited(const QString &arg1);
    void on_initEnergy_textEdited(const QString &arg1);
    void on_targ_element_textEdited(const QString &arg1);
    void on_targ_A_textEdited(const QString &arg1);

    void on_actionExecute_triggered();
    void on_actionAbout_triggered();
    void on_action_Exit_triggered();
    void on_actionPrint_triggered();
    void on_actionWeb_Documentation_triggered();

    void on_actionOpen_triggered();
    void on_actionSave_triggered();
    void on_actionSave_As_triggered();

    void on_actionPrint_Preview_triggered();
    void printPreview(QPrinter*);

    void on_actionNIM_triggered();

private:
    Ui::MainWindow *ui;

};

#endif // MAINWINDOW_H
//------------------------------------------------
