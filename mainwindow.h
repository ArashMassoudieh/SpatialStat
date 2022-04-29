#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "environment.h"

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
    Environment environment;
private:
    Ui::MainWindow *ui;

public slots:
    void on_test();
};
#endif // MAINWINDOW_H
