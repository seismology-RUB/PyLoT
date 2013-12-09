#include "qt_pylot.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    Qt_PyLoT w;
    w.show();

    return a.exec();
}
