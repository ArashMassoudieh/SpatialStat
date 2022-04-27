#include "command.h"
#include "vector"
#include "Utilities.h"

vector<char> Command::deliminators = vector<char>({'=', '(', ')'});


Command::Command()
{

}

Command::Command(const string &cmdline)
{
    if (aquiutils::contains(cmdline,"."))
    {
        CommadType = command_type::modifier;
    }
    else if (aquiutils::contains(cmdline,"="))
    {
        CommadType = command_type::creator;
    }
    else
        return;
    vector<string> s = aquiutils::split(cmdline,deliminators);
    if (s.size()>=2)
    {
        command = s[1];
        object_name = s[0];
        for (unsigned int i=2; i<s.size(); i++)
        {
            vector<string> arg = aquiutils::split(s[i],':');
            if (arg.size()==2)
                arguments[arg[0]] = arg[1];
        }
    }
}

