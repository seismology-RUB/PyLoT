#include "qt_pylot.h"
#include "ui_qt_pylot.h"

Qt_PyLoT::Qt_PyLoT(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::Qt_PyLoT)
{
    ui->setupUi(this);
}

Qt_PyLoT::~Qt_PyLoT()
{
    delete ui;
}
