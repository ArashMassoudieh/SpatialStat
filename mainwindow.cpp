#include "mainwindow.h"
#include "ui_mainwindow.h"


MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    connect(ui->actionRun_Sample_Script,SIGNAL(triggered()),this,SLOT(on_test()));
    ui->tableWidget->setColumnCount(3);
    ui->tableWidget->setRowCount(0);
    ui->tableWidget->setHorizontalHeaderItem(0,new QTableWidgetItem(QString("Command")));
    ui->tableWidget->setHorizontalHeaderItem(1,new QTableWidgetItem(QString("Parameters")));
    ui->tableWidget->setHorizontalHeaderItem(2,new QTableWidgetItem(QString("Progress")));
    environment.outputwindow = ui->tableWidget;
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_test()
{
    Command cmd = Command("grid=CreateGrid(nx:100,ny:100,dx:0.01,dy:0.01)");
    environment.Execute(cmd);
    cmd = Command("dist=CreateDistribution(type:lognormal,p0:1,p1:0,p2:1)");
    environment.Execute(cmd);
    cmd = Command("dist*WriteToFile(filename:test.txt,nbins:50");
    environment.Execute(cmd);
    cmd = Command("dist*SetInverseCumulative(ninc:1000)");
    environment.Execute(cmd);
    cmd = Command("dist*WriteInverseCumulativeToFile(filename:inversecumulative.txt)");
    environment.Execute(cmd);
    cmd = Command("grid*AssignKField(Distribution:dist,correlation_length_x:0.5,correlation_length_y:0.5,Maximum_neighboring_nodes:11)");
    environment.Execute(cmd);
    cmd = Command("grid*RenormalizeKField(Distribution:dist");
    environment.Execute(cmd);
    cmd = Command("grid*WriteKFieldToVTP(filename:K_field.vtp,z_scale:0.1,log_scale:0");
    environment.Execute(cmd);
    cmd = Command("grid*SolveHydro(l_boundary:1,r_boundary=0");
    environment.Execute(cmd);
    cmd = Command("grid*WriteHydroSolutionToVTP(filename:hydro_solution.vtp,z_scale:0.1,log_scale:0");
    environment.Execute(cmd);
    cmd = Command("grid*SolveTransport(nspecies:1,decay_coeff:1,decay_order:1,time_weight:1,l_boundary:1,diffusion:0,dt:0.005,t_end:1)");
    environment.Execute(cmd);


}
