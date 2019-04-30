function Table(n) {
    this.length=n;
    return this; 
	}

function DateModif() {
    NameMonth   =new Table(12);
    NameMonth[1] ="January";
    NameMonth[2] ="February";
    NameMonth[3] ="Mars";
    NameMonth[4] ="April";
    NameMonth[5] ="Mai";
    NameMonth[6] ="June";
    NameMonth[7] ="Jully"; 
    NameMonth[8] ="August";
    NameMonth[9] ="September";
    NameMonth[10]="October";
    NameMonth[11]="November";
    NameMonth[12]="December";
    Date       =new Date(document.lastModified);
    var Month   =NameMonth[Date.getMonth()+1];
    var Year  =Date.getFullYear();
    return Date.getDate()+" "+Month+" "+Year;
	}