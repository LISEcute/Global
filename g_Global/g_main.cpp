#include "g_mainwindow.h"
#include <QApplication>
#include <QStyleFactory>

#include <QDir>
#include <QSettings>
#include <QStandardPaths>

extern void getInitialDir(void);
extern void getDPIscaling(void);

extern QString LISErootPATH;
extern QString MyDocCompPATH;
extern const char *LISEini;

QString FileArg="";
int fontsizeGlobal=9;
int useHighDpiScaling=1;

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

int main(int argc, char *argv[])
{
//getDPIscaling();

//if(!useHighDpiScaling)
//  QCoreApplication::setAttribute(Qt::AA_Use96Dpi);

//   QApplication::setAttribute(useHighDpiScaling ? Qt::AA_EnableHighDpiScaling :
//                                                   Qt::AA_DisableHighDpiScaling);
//--------------------------------------------------------------------------------------
    QApplication app(argc, argv);
    LISErootPATH = QCoreApplication::applicationDirPath();
    if(LISErootPATH.size()<=0)
        {
        QString arg0 = QDir::fromNativeSeparators(argv[0]);
        QDir dir(arg0);
        LISErootPATH = dir.currentPath();
        }

    if(argc>1) FileArg = QDir::fromNativeSeparators(argv[1]);

    //---------------------------------------------------------------------------------- paths + read font size start
    getInitialDir();
    //---------------------------------------------------------------------------------- paths + read font size stop

    //------------------------------------------------------------------ style start

    app.setFont(QFont("Arial", fontsizeGlobal, QFont::Normal));
    app.setOrganizationName("FRIB/MSU");
    app.setStyle(QStyleFactory::create("Fusion"));
    QFile FileStyle(":/w_Main/mainstyle.qss");
    FileStyle.open(QFile::ReadOnly);
    QString StyleSheet = QLatin1String(FileStyle.readAll());
    qApp->setStyleSheet(StyleSheet);
     //------------------------------------------------------------------ style end

    MainWindow w;
    w.show();

    return app.exec();
}
