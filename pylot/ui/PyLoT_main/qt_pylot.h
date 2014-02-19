#ifndef QT_PYLOT_H
#define QT_PYLOT_H

#include <QMainWindow>

namespace Ui {
class Qt_PyLoT;
}

class Qt_PyLoT : public QMainWindow
{
    Q_OBJECT

public:
    explicit Qt_PyLoT(QWidget *parent = 0);
    ~Qt_PyLoT();

private:
    Ui::Qt_PyLoT *ui;
};

#endif // QT_PYLOT_H
