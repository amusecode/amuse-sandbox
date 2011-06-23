import pygtk
pygtk.require('2.0')
import gtk
import pango
import mathmod

numpad_buttons=( (0,1,3,4,'0'),
                 (1,2,3,4,'.'),
                 (0,1,2,3,'1'),
                 (1,2,2,3,'2'),
                 (2,3,2,3,'3'),
                 (0,1,1,2,'4'),
                 (1,2,1,2,'5'),
                 (2,3,1,2,'6'),
                 (0,1,0,1,'7'),
                 (1,2,0,1,'8'),
                 (2,3,0,1,'9'),
                 (3,4,0,1,'+'),
                 (3,4,1,2,'-'),
                 (3,4,2,3,'*'),
                 (3,4,3,4,'/'),
                 (4,5,1,2,'('),
                 (5,6,1,2,')') )

shiftable_buttons=((1,2,0,1,'sin','sin(','asin','asin('),
                   (2,3,0,1,'cos','cos(','acos','acos('),
                   (3,4,0,1,'tan','tan(','atan','atan('),
                   (1,2,1,2,'sinh','sinh(','asinh','asinh('),
                   (2,3,1,2,'cosh','cosh(','acosh','acosh('),
                   (3,4,1,2,'tanh','tanh(','atanh','atanh('))

scipad_buttons=((2,3,2,3,'x^2','**2'),
                (3,4,2,3,'sqrt','sqrt('),
                (4,5,2,3,'x^-1','**-1'),
                (6,7,0,1,'|','|'),
                (6,7,1,2,'abs','abs('),
                (6,7,2,3,'rnd','rnd()'),
                (0,1,1,2,'pi','pi'),
                (4,5,0,1,'log','log('),
                (4,5,1,2,'log10','log10('),                
                (5,6,0,1,'e','e'),
                (5,6,1,2,'10^x','10**('),
                (5,6,2,3,'^','**'))

unitpad_buttons=((0,1,1,2,'convert',').convert_to('),
                 (1,2,1,2,'value',').value_in('),
                 (2,3,1,2,'in base','in_base('),
                 (3,4,1,2,'strip','strip('),
                 (4,5,1,2,'clean','clean('),
                 (5,6,1,2,'simplify','simplify('))

special_functions=('atan2','fac','cot','coth')

si_prefixes=mathmod.si_prefixes()
constants=mathmod.constants()
siunits=mathmod.siunits()
named_units=mathmod.named_units()
                 
class unitcalc(object):
  def insert_in_buff(self,data):
    self.buff.insert_at_cursor(data)
    self.view.scroll_mark_onscreen(self.buff.get_insert())

  def token_callback(self,widget,data):
    self.insert_in_buff(data)

  def bksp_callback(self,widget):
    start=self.buff.get_iter_at_line(-1)
    end=self.buff.get_end_iter()
    current=self.buff.get_iter_at_mark(self.buff.get_insert())
    if start.compare(current):
      self.buff.backspace(current,True,True)
    
  def clr_callback(self,widget):
    start=self.buff.get_iter_at_line(-1)
    end=self.buff.get_end_iter()
    current=self.buff.get_iter_at_mark(self.buff.get_insert())
    if start.compare(current)==0 and end.compare(current)==0:
      start=self.buff.get_iter_at_line(0)
    self.buff.delete(start,end)


  def evaluate_callback(self,widget):
    start=self.buff.get_iter_at_line(-1)
    end=self.buff.get_end_iter()
    express=self.buff.get_text(start,end)
    unclosedleft=0
    unclosedright=0
    for x in express:
      if(x==')'):
        if(unclosedleft==0):
          unclosedright+=1
        else:
          unclosedleft+=-1  
      if(x=='('):
        unclosedleft+=1  
    self.buff.insert(end, unclosedleft*')')
    start=self.buff.get_iter_at_line(-1)
    end=self.buff.get_end_iter()
    self.buff.insert(start, unclosedright*'(')
    start=self.buff.get_iter_at_line(-1)
    end=self.buff.get_end_iter()
    express=self.buff.get_text(start,end)
    self.buff.insert(end,'=\n')
    end=self.buff.get_end_iter()
    try:
      result=self.calcengine.calculate(express)
    except ZeroDivisionError:
      self.buff.insert_with_tags_by_name(end,"math error\n","right_just")
      end=self.buff.get_end_iter()
      self.buff.insert(end,express)
      self.view.scroll_mark_onscreen(self.buff.get_insert())      
    except (SyntaxError,TypeError):
      self.buff.insert_with_tags_by_name(end,"syntax error\n","right_just")
      end=self.buff.get_end_iter()
      self.buff.insert(end,express)
      self.view.scroll_mark_onscreen(self.buff.get_insert())      
#    except:
#      self.buff.insert_with_tags_by_name(end,"error\n","right_just")
#      end=self.buff.get_end_iter()
#      self.buff.insert(end,express)
#      self.view.scroll_mark_onscreen(self.buff.get_insert())
    else:
      self.calcengine.pushback_result(result)
      self.calcengine.pushback_express(express)      
      self.buff.insert_with_tags_by_name(end,str(result)+"\n","right_just")
      self.view.scroll_mark_onscreen(self.buff.get_insert())

    

  def move_callback(self,widget,data):
    start=self.buff.get_iter_at_line(-1)
    end=self.buff.get_end_iter()
    current=self.buff.get_iter_at_mark(self.buff.get_insert())
    if(data=='<'):
#      if start.compare(current):
#        current.backward_cursor_position() 
      self.buff.delete(start,end)
      express=self.calcengine.backward_express()
      self.insert_in_buff(express)
      current=self.buff.get_iter_at_mark(self.buff.get_insert())
    if(data=='>'):
#      if current.compare(end):
#        current.forward_cursor_position()
      self.buff.delete(start,end)
      express=self.calcengine.forward_express()
      self.insert_in_buff(express)
      current=self.buff.get_iter_at_mark(self.buff.get_insert())
    self.buff.place_cursor(current)



  def ans_callback(self,widget):
    self.insert_in_buff("ans")
    

  def shift_callback(self,widget):
    self.shift=not self.shift
    if(self.shift):
      l=6
    else:
      l=4
    for i,x in enumerate(self.trigobuttons):
      x.set_label(shiftable_buttons[i][l])

  def trigo_callback(self,widget,data):
    if(self.shift):
      self.insert_in_buff('a'+data)
      self.shift_button.set_active(False)
    else:
      self.insert_in_buff(data)

  def combobox_callback(self,widget,data):
    model=widget.get_model()
    index=widget.get_active()
    if(widget.ignore.count(index) == 0):
      self.insert_in_buff(model[index][0]+data)      
      widget.set_active(0)
      
                
  def __init__(self):
    self.window = gtk.Window(gtk.WINDOW_TOPLEVEL)

    self.window.set_title("AMUSE Calculator")
    self.window.set_size_request(480,640)
  
    self.window.connect("delete_event", gtk.main_quit)

    self.table = gtk.Table(12,1,True) 
    self.window.add(self.table)
    
    self.buff=gtk.TextBuffer()
    self.view=gtk.TextView(self.buff)
    self.view.set_editable(False)
    self.view.set_cursor_visible(True)
    self.view.set_wrap_mode(gtk.WRAP_CHAR)
    self.view.modify_base(gtk.STATE_NORMAL, gtk.gdk.color_parse("dark gray"))
    self.view.modify_font(pango.FontDescription("Monospace 14"))
    self.view.set_border_width(3)
    self.view.set_left_margin(3)
    self.view.set_right_margin(3)
    self.table.attach(self.view,0,1,0,3)

    self.numpad=gtk.Table(4,6,True)
    self.scipad=gtk.Table(3,7,True)   
    self.unitpad=gtk.Table(2,6,True)

    self.table.attach(self.unitpad,0,1,3,5)
    self.table.attach(self.scipad,0,1,5,8)
    self.table.attach(self.numpad,0,1,8,12)
    
    for x in numpad_buttons:
      button=gtk.Button(x[4])
      button.connect("clicked",self.token_callback,x[4])
      self.numpad.attach(button,x[0],x[1],x[2],x[3])
    
    button=gtk.Button('off')
    button.connect("clicked",gtk.main_quit)
    self.numpad.attach(button,5,6,0,1)

    button=gtk.Button('=')
    button.connect("clicked",self.evaluate_callback)
    self.numpad.attach(button,2,3,3,4)

    button=gtk.Button('bksp')
    button.connect("clicked",self.bksp_callback)
    self.numpad.attach(button,4,5,3,4)

    button=gtk.Button('clr')
    button.connect("clicked",self.clr_callback)
    self.numpad.attach(button,5,6,3,4)

    button=gtk.Button('<')
    button.connect("clicked",self.move_callback,'<')
    self.numpad.attach(button,4,5,2,3)
    button=gtk.Button('>')
    button.connect("clicked",self.move_callback,'>')
    self.numpad.attach(button,5,6,2,3)

    button=gtk.Button('ans')
    button.connect("clicked",self.ans_callback)
    self.numpad.attach(button,4,5,0,1)

    for x in scipad_buttons:
      button=gtk.Button(x[4])
      button.connect("clicked",self.token_callback,x[5])
      self.scipad.attach(button,x[0],x[1],x[2],x[3])

    self.func_combo=gtk.combo_box_new_text()
    self.func_combo.ignore=[0]
    self.func_combo.append_text("special func")
    for x in special_functions:
      self.func_combo.append_text(x)
    self.scipad.attach(self.func_combo,0,2,2,3)
    self.func_combo.connect('changed',self.combobox_callback,'(')
    self.func_combo.set_active(0)


    self.shift_button=gtk.ToggleButton('shift')
    self.shift_button.connect("clicked",self.shift_callback)
    self.scipad.attach(self.shift_button,0,1,0,1)

    self.trigobuttons=[]
    for x in shiftable_buttons:
      button=gtk.Button(x[4])
      button.connect("clicked",self.trigo_callback,x[5])
      self.scipad.attach(button,x[0],x[1],x[2],x[3])
      self.trigobuttons.append(button)
    
    for x in unitpad_buttons:
      button=gtk.Button(x[4])
      button.connect("clicked",self.token_callback,x[5])
      self.unitpad.attach(button,x[0],x[1],x[2],x[3])

    self.constants_combo=gtk.combo_box_new_text()
    self.constants_combo.ignore=[0]
    self.constants_combo.append_text("constants")
    for x in constants:
      self.constants_combo.append_text(x)
    self.unitpad.attach(self.constants_combo,4,6,0,1)
    self.constants_combo.connect('changed',self.combobox_callback,'')
    self.constants_combo.set_active(0)

    self.siunits_combo=gtk.combo_box_new_text()
    self.siunits_combo.ignore=[0,8]
    self.siunits_combo.append_text("SI units")
    for x in siunits:
      self.siunits_combo.append_text(x)
    self.siunits_combo.append_text("--- prefixes ---")
    for x in si_prefixes:
      self.siunits_combo.append_text(x)
    self.unitpad.attach(self.siunits_combo,0,2,0,1)
    self.siunits_combo.connect('changed',self.combobox_callback,'')
    self.siunits_combo.set_active(0)

    self.named_combo=gtk.combo_box_new_text()
    self.named_combo.ignore=[0]
    self.named_combo.append_text("misc units")
    for x in named_units:
      self.named_combo.append_text(x)
    self.unitpad.attach(self.named_combo,2,4,0,1)
    self.named_combo.connect('changed',self.combobox_callback,'')
    self.named_combo.set_active(0)


    self.window.show_all()

    self.buff.create_tag("right_just", justification=gtk.JUSTIFY_RIGHT)
    self.shift=False

    self.calcengine=mathmod.calcengine()
        
def main():
    gtk.main()

if __name__ == "__main__":
    hello = unitcalc()
    main()
    
    
    
