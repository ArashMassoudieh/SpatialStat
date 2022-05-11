#include "interface.h"
#include "Utilities.h"
#include "environment.h"
#include "QCoreApplication"

Interface::Interface(Environment *_parent)
{
    parent = _parent;
}

bool Interface::HasCommand(const string &cmd)
{
    if (aquiutils::lookup(commands(),cmd)!=-1)
        return true;
    else
        return false;
}

vector<string> Interface::commands()
{
    return vector<string>();
}

QTableWidget *Interface::outputwindow()
{
    return parent->outputwindow;
}

void Interface::SetProgressValue(const double &x)
{
    parent->current_progress_bar->setValue(x*100);
    outputwindow()->update();
    QCoreApplication::processEvents();
}
