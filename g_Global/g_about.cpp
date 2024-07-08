#include "g_about.h"
#include "ui_g_about.h"
#include "g_ftype.h"
#include <QDesktopServices>
#include <QUrl>

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
About::About(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::About)
{
    ui->setupUi(this);
    ui->label_Version->setText(Global_version);
    ui->label_Date->setText(Global_date);

    connect(ui->label_LISE, SIGNAL(clicked()), this, SLOT(CmLISE()));
    connect(ui->label_Charge, SIGNAL(clicked()), this, SLOT(CmCharge()));
    connect(ui->label_mail, SIGNAL(clicked()), this, SLOT(CmMail()));

}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
About::~About() {    delete ui; }
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void About::CmLISE() {    QDesktopServices::openUrl(QUrl(ui->label_LISE->text()));}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void About::CmCharge() {    QDesktopServices::openUrl(QUrl(ui->label_Charge->text()));}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void About::CmMail()
{
QString ss("mailto:tarasov@frib.msu.edu?subject=Global ");
ss.append(ui->label_Version->text());
QDesktopServices::openUrl(QUrl(ss));
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
