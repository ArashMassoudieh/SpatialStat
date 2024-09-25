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
    Command cmd1 = Command("grid=CreateGrid(nx=100,ny=100,dx=0.01,dy=0.01)");
    environment.Execute(cmd1);
    Command cmd2 = Command("dist=CreateDistribution(type=lognormal,p0=1,p1=0,p2=1)");
    environment.Execute(cmd2);
    Command cmd3 = Command("dist.WriteDistToFile(filename=test.txt,nbins=50");
    environment.Execute(cmd3);
    Command cmd4 = Command("dist.SetInverseCumulative(ninc=1000)");
    environment.Execute(cmd4);
    Command cmd5 = Command("dist.WriteInverseCumulativeToFile(filename=inversecumulative.txt)");
    environment.Execute(cmd5);
    Command cmd6 = Command("grid.AssignKField(Distribution=dist,correlation_length_x=0.5,correlation_length_y=0.5,Maximum_neighboring_nodes=11)");
    environment.Execute(cmd6);
    Command cmd7 = Command("grid.RenormalizeKField(Distribution=dist");
    environment.Execute(cmd7);
    Command cmd8 = Command("grid.WriteKFieldToVTP(filename=K_field.vtp,z_scale=0.1,log_scale=0");
    environment.Execute(cmd8);
    Command cmd9 = Command("grid.SolveHydro(l_boundary=1,r_boundary=0");
    environment.Execute(cmd9);
    Command cmd10 = Command("grid.WriteHydroSolutionToVTP(filename=hydro_solution.vtp,z_scale=0.1,log_scale=0");
    environment.Execute(cmd10);
    Command cmd11 = Command("grid.SolveTransport(nspecies=1,decay_coeff=1,decay_order=1,time_weight=1,l_boundary=1,diffusion=0,dt=0.005,t_end=1)");
    environment.Execute(cmd11);
    Command cmd12 = Command("grid.WriteConcentrationToVTP(filename=Concentration.vtp,interval=5)");
    environment.Execute(cmd12);
    Command cmd13 = Command("BTC1=grid.GetConcentrationBTCAtX(x=0.5,filename=BTC_c.txt,filename_d=BTC.txt)");
    environment.Execute(cmd13);
    Command cmd14 = Command("pthways1=grid.CreateTrajectories(n=20,x_0=0.01,dx=0.01,x_end=0.95,tol=0.001)");
    environment.Execute(cmd14);
    Command cmd15 = Command("pthways1.Uniformize(dx=0.01)");
    environment.Execute(cmd15);
    Command cmd16 = Command("pthways1.WritePathwayToVTP(filename=paths.vtp)");
    environment.Execute(cmd16);
    Command cmd17 = Command("vdist = grid.GetMarginalVelocityDistribution(direction=x)");
    environment.Execute(cmd17);
    Command cmd18 = Command("vdist.WriteTimeSeriesToFile(filename=v_dist.txt,nbins=50");
    environment.Execute(cmd18);
    Command cmd19 = Command("grid_2nd=grid.AssignKField_2nd_order(correlation_length_x=10,correlation_length_y=0.1)");
    environment.Execute(cmd19);
    Command cmd20 = Command("grid_2nd.WriteFieldToVTP(filename=2nd_order_k_score.vtp,prop=K_Gauss");
    environment.Execute(cmd20);
    Command cmd21 = Command("dist_K2=grid_2nd.GetMarginalDistribution(prop=K_Gauss, nbins=30");
    environment.Execute(cmd21);
    Command cmd22 = Command("dist_K2.WriteDistToFile(filename=distK2.txt,nbins=30");
    environment.Execute(cmd22);



}
