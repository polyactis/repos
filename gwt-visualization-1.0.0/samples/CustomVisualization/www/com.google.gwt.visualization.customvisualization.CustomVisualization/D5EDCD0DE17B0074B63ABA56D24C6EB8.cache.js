(function(){var $gwt_version = "1.5.3";var $wnd = window;var $doc = $wnd.document;var $moduleName, $moduleBase;var $stats = $wnd.__gwtStatsEvent ? function(a) {return $wnd.__gwtStatsEvent(a);} : null;$stats && $stats({moduleName:'com.google.gwt.visualization.customvisualization.CustomVisualization',subSystem:'startup',evtGroup:'moduleStartup',millis:(new Date()).getTime(),type:'moduleEvalStart'});var B='',F=', Row size: ',mb=', Size: ',ub=':',sb='<\/b>',rb='<b>',gb='Cannot create a column with a negative index: ',hb='Cannot create a row with a negative index: ',Ab='CustomVisualization',v='DOMMouseScroll',lb='Index: ',Eb='Integer;',bc='JavaScriptObject$;',z='MouseEvents',Fb='Object;',E='Row index: ',Cb='Widget;',ac='[Lcom.google.gwt.core.client.',Bb='[Lcom.google.gwt.user.client.ui.',Db='[Ljava.lang.',kb='__widgetID',pb='backgroundColor',l='blur',m='change',x='click',jb='col',ib='colgroup',zb='com.google.gwt.visualization.customvisualization.client.CustomVisualizationEntryPoint',w='contextmenu',cb='dblclick',t='error',nb='focus',y='html',yb='keydown',cc='keypress',dc='keyup',A='left',ec='load',fc='losecapture',wb='moduleStartup',n='mousedown',o='mousemove',p='mouseout',q='mouseover',r='mouseup',u='mousewheel',xb='onModuleLoadStart',D='position',s='scroll',ob='select',tb='selection changed',vb='startup',eb='table',ab='tagName',fb='tbody',bb='td',C='top',db='tr',qb='white';var _;function an(a){return (this==null?null:this)===(a==null?null:a)}
function bn(){return this.$H||(this.$H=++yc)}
function Em(){}
_=Em.prototype={};_.eQ=an;_.hC=bn;_.tM=as;_.tI=1;function mc(b,a){return b.tM==as||b.tI==2?b.eQ(a):(b==null?null:b)===(a==null?null:a)}
function oc(a){return a.tM==as||a.tI==2?a.hC():a.$H||(a.$H=++yc)}
var yc=0;function Cc(b){var a=b.firstChild;while(a&&a.nodeType!=1)a=a.nextSibling;return a}
function Dc(a){var b=a.parentNode;if(b==null){return null}if(b.nodeType!=1)b=null;return b}
function Ec(a,b){while(a.firstChild){a.removeChild(a.firstChild)}if(b!=null){a.appendChild($doc.createTextNode(b))}}
function ud(e,c){var d=[null,0,false,[0,0]];var f=d[e];var a=new Array(c);for(var b=0;b<c;++b){a[b]=f}return a}
function vd(a,f,c,b,e){var d;d=ud(e,b);wd(a,f,c,d);return d}
function wd(b,d,c,a){if(!xd){xd=new qd()}Ad(a,xd);a.tI=d;a.qI=c;return a}
function yd(a,b,c){if(c!=null){if(a.qI>0&&!Dd(c.tI,a.qI)){throw new Dl()}if(a.qI<0&&(c.tM==as||c.tI==2)){throw new Dl()}}return a[b]=c}
function Ad(a,c){for(var b in c){var d=c[b];if(d){a[b]=d}}return a}
function qd(){}
_=qd.prototype=new Em();_.tI=0;_.length=0;_.qI=0;var xd=null;function Ed(b,a){return b&&!!ke[b][a]}
function Dd(b,a){return b&&ke[b][a]}
function Fd(b,a){if(b!=null&&!Dd(b.tI,a)){throw new bm()}return b}
function ce(b,a){return b!=null&&Ed(b.tI,a)}
var ke=[{},{},{1:1,17:1,18:1,19:1},{2:1,17:1},{2:1,17:1},{2:1,17:1},{2:1,17:1},{9:1},{4:1,5:1,9:1},{4:1,5:1,8:1,9:1},{4:1,5:1,8:1,9:1},{4:1,5:1,8:1,9:1},{4:1,5:1,9:1},{4:1,5:1,8:1,9:1},{4:1,5:1,8:1,9:1},{4:1,5:1,6:1,8:1,9:1},{3:1},{4:1,5:1,6:1,8:1,9:1},{13:1},{13:1,17:1},{13:1,17:1},{4:1,5:1,9:1},{4:1,5:1,9:1},{7:1},{2:1,17:1},{2:1,17:1},{2:1,17:1},{2:1,17:1},{2:1,17:1},{17:1,20:1},{2:1,17:1},{2:1,17:1},{11:1,17:1,19:1,20:1},{18:1},{18:1},{2:1,17:1},{15:1},{15:1},{12:1},{12:1},{12:1},{15:1},{14:1,17:1},{15:1,17:1},{12:1},{2:1,17:1},{10:1}];function bf(b,a,c){var d;if(a==ef){if(hg(b)==8192){ef=null}}d=af;af=b;try{c.z(b)}finally{af=d}}
function ff(a,b){jg();bg(a,b);Af(a,b)}
var af=null,ef=null;function lf(a){sf();if(!nf){nf=uq(new tq())}vq(nf,a)}
function of(){var a;if(nf){for(a=lp(new jp(),nf);a.a<a.b.ab();){Fd(op(a),3);hj()}}}
function pf(){var a,b;b=null;if(nf){for(a=lp(new jp(),nf);a.a<a.b.ab();){Fd(op(a),3);b=null}}return b}
function rf(){__gwt_initHandlers(function(){},function(){return pf()},function(){of()})}
function sf(){if(!qf){rf();qf=true}}
var nf=null,qf=false;function hg(a){switch(a.type){case l:return 4096;case m:return 1024;case x:return 1;case cb:return 2;case nb:return 2048;case yb:return 128;case cc:return 256;case dc:return 512;case ec:return 32768;case fc:return 8192;case n:return 4;case o:return 64;case p:return 32;case q:return 16;case r:return 8;case s:return 16384;case t:return 65536;case u:return 131072;case v:return 131072;case w:return 262144;}}
function jg(){if(!lg){Ff();zf();lg=true}}
function mg(a){return a!=null&&Ed(a.tI,4)&&!(a!=null&&(a.tM!=as&&a.tI!=2))}
var lg=false;function Ef(c,e){var b=0,a=c.firstChild;while(a){if(a===e){return b}if(a.nodeType==1){++b}a=a.nextSibling}return -1}
function Ff(){eg=function(b){if(dg(b)){var a=cg;if(a&&a.__listener){if(mg(a.__listener)){bf(b,a,a.__listener);b.stopPropagation()}}}};dg=function(a){return true};fg=function(b){var c,a=this;while(a&&!(c=a.__listener)){a=a.parentNode}if(a&&a.nodeType!=1){a=null}if(c){if(mg(c)){bf(b,a,c)}}};$wnd.addEventListener(x,eg,true);$wnd.addEventListener(cb,eg,true);$wnd.addEventListener(n,eg,true);$wnd.addEventListener(r,eg,true);$wnd.addEventListener(o,eg,true);$wnd.addEventListener(q,eg,true);$wnd.addEventListener(p,eg,true);$wnd.addEventListener(u,eg,true);$wnd.addEventListener(yb,dg,true);$wnd.addEventListener(dc,dg,true);$wnd.addEventListener(cc,dg,true)}
function ag(e,g,d){var c=0,b=e.firstChild,a=null;while(b){if(b.nodeType==1){if(c==d){a=b;break}++c}b=b.nextSibling}e.insertBefore(g,a)}
function bg(c,a){var b=(c.__eventBits||0)^a;c.__eventBits=a;if(!b)return;if(b&1)c.onclick=a&1?fg:null;if(b&2)c.ondblclick=a&2?fg:null;if(b&4)c.onmousedown=a&4?fg:null;if(b&8)c.onmouseup=a&8?fg:null;if(b&16)c.onmouseover=a&16?fg:null;if(b&32)c.onmouseout=a&32?fg:null;if(b&64)c.onmousemove=a&64?fg:null;if(b&128)c.onkeydown=a&128?fg:null;if(b&256)c.onkeypress=a&256?fg:null;if(b&512)c.onkeyup=a&512?fg:null;if(b&1024)c.onchange=a&1024?fg:null;if(b&2048)c.onfocus=a&2048?fg:null;if(b&4096)c.onblur=a&4096?fg:null;if(b&8192)c.onlosecapture=a&8192?fg:null;if(b&16384)c.onscroll=a&16384?fg:null;if(b&32768)c.onload=a&32768?fg:null;if(b&65536)c.onerror=a&65536?fg:null;if(b&131072)c.onmousewheel=a&131072?fg:null;if(b&262144)c.oncontextmenu=a&262144?fg:null}
var cg=null,dg=null,eg=null,fg=null;function zf(){$wnd.addEventListener(p,function(b){var a=$wnd.__captureElem;if(a&&!b.relatedTarget){if(y==b.target.tagName.toLowerCase()){var c=$doc.createEvent(z);c.initMouseEvent(r,true,true,$wnd,0,b.screenX,b.screenY,b.clientX,b.clientY,b.ctrlKey,b.altKey,b.shiftKey,b.metaKey,b.button,null);a.dispatchEvent(c)}}},true);$wnd.addEventListener(v,eg,true)}
function Af(b,a){if(a&131072){b.addEventListener(v,fg,false)}}
function uj(b,a){b.i=a}
function sj(){}
_=sj.prototype=new Em();_.tI=7;_.i=null;function fk(a){if(a.v()){throw new mm()}a.g=true;a.i.__listener=a;a.m();a.B()}
function gk(a){if(!a.v()){throw new mm()}try{a.C()}finally{a.n();a.i.__listener=null;a.g=false}}
function hk(a){if(ce(a.h,8)){Fd(a.h,8).D(a)}else if(a.h){throw new mm()}}
function ik(c,b){var a;a=c.h;if(!b){if(!!a&&a.v()){c.A()}c.h=null}else{if(a){throw new mm()}c.h=b;if(b.v()){c.y()}}}
function jk(){}
function kk(){}
function lk(){return this.g}
function mk(){fk(this)}
function nk(a){}
function ok(){gk(this)}
function pk(){}
function qk(){}
function vj(){}
_=vj.prototype=new sj();_.m=jk;_.n=kk;_.v=lk;_.y=mk;_.z=nk;_.A=ok;_.B=pk;_.C=qk;_.tI=8;_.g=false;_.h=null;function Bi(){var a,b;for(b=this.w();b.u();){a=Fd(b.x(),5);a.y()}}
function Ci(){var a,b;for(b=this.w();b.u();){a=Fd(b.x(),5);a.A()}}
function Di(){}
function Ei(){}
function zi(){}
_=zi.prototype=new vj();_.m=Bi;_.n=Ci;_.B=Di;_.C=Ei;_.tI=9;function ug(c,a,b){hk(a);Ej(c.a,a);b.appendChild(a.i);ik(a,c)}
function wg(b,c){var a;if(c.h!=b){return false}ik(c,null);a=c.i;Dc(a).removeChild(a);dk(b.a,c);return true}
function xg(){return zj(new xj(),this.a)}
function yg(a){return wg(this,a)}
function sg(){}
_=sg.prototype=new zi();_.w=xg;_.D=yg;_.tI=10;function og(a,b){ug(a,b,a.i)}
function qg(a){a.style[A]=B;a.style[C]=B;a.style[D]=B}
function rg(b){var a;a=wg(this,b);if(a){qg(b.i)}return a}
function ng(){}
_=ng.prototype=new sg();_.D=rg;_.tI=11;function Bg(a,b){if(a.e){throw new mm()}hk(b);uj(a,b.i);a.e=b;ik(b,a)}
function Cg(){if(this.e){return this.e.g}return false}
function Dg(){fk(this.e);this.i.__listener=this}
function Eg(a){oi(this.e,a)}
function Fg(){gk(this.e)}
function zg(){}
_=zg.prototype=new vj();_.v=Cg;_.y=Dg;_.z=Eg;_.A=Fg;_.tI=12;_.e=null;function di(b,a){if(!b.e){b.e=oj(new nj());ff(b.i,1|(b.i.__eventBits||0))}vq(b.e,a)}
function ei(c,a){var b;b=c.a.rows.length;if(a>=b||a<0){throw qm(new pm(),E+a+F+b)}}
function gi(d){var a,b,c;for(c=0;c<d.a.rows.length;++c){for(b=0;b<(ei(d,c),d.a.rows[c].cells.length);++b){a=li(d,c,b);if(a){pi(d,a)}}}}
function ki(d,b){var a,c,e;c=b.target;for(;c;c=Dc(c)){if(qn(c[ab]==null?null:String(c[ab]),bb)){e=Dc(c);a=Dc(e);if(a==d.a){return c}}if(c==d.a){return null}}return null}
function li(e,d,b){var a,c;c=e.b.a.a.rows[d].cells[b];a=Cc(c);if(!a){return null}else{return Dh(e.f,a)}}
function mi(b,a){var c;if(a!=b.a.rows.length){ei(b,a)}c=$doc.createElement(db);ag(b.a,c,a);return a}
function ni(d,c,a){var b,e;b=Cc(c);e=null;if(b){e=Dh(d.f,b)}if(e){pi(d,e);return true}else{if(a){c.innerHTML=B}return false}}
function oi(f,c){var a,b,d,e,g;switch(hg(c)){case 1:{if(f.e){e=ki(f,c);if(!e){return}g=Dc(e);a=Dc(g);d=Ef(a,g);b=Ef(g,e);qj(f.e,d,b)}break}}}
function pi(b,c){var a;if(c.h!=b){return false}ik(c,null);a=c.i;Dc(a).removeChild(a);Eh(b.f,a);return true}
function ri(b,a){b.c=a;rh(b.c)}
function si(f,d,a,c){var e,b;gh(f,d,a);e=(b=f.b.a.a.rows[d].cells[a],ni(f,b,c==null),b);if(c!=null){e.innerHTML=c||B}}
function ti(f,c,a,e){var d,b;gh(f,c,a);d=(b=f.b.a.a.rows[c].cells[a],ni(f,b,e==null),b);if(e!=null){Ec(d,e)}}
function ui(){return vh(new th(),this.f)}
function vi(a){oi(this,a)}
function wi(a){return pi(this,a)}
function jh(){}
_=jh.prototype=new zi();_.w=ui;_.z=vi;_.D=wi;_.tI=13;_.a=null;_.b=null;_.c=null;_.d=null;_.e=null;function eh(a){a.f=Bh(new sh());a.d=$doc.createElement(eb);a.a=$doc.createElement(fb);a.d.appendChild(a.a);a.i=a.d;a.b=ch(new bh(),a);ri(a,ph(new oh(),a));return a}
function gh(e,d,b){var a,c;hh(e,d);if(b<0){throw qm(new pm(),gb+b)}a=(ei(e,d),e.a.rows[d].cells.length);c=b+1-a;if(c>0){ih(e.a,d,c)}}
function hh(d,b){var a,c;if(b<0){throw qm(new pm(),hb+b)}c=d.a.rows.length;for(a=c;a<=b;++a){mi(d,a)}}
function ih(f,d,c){var e=f.rows[d];for(var b=0;b<c;b++){var a=$doc.createElement(bb);e.appendChild(a)}}
function ah(){}
_=ah.prototype=new jh();_.tI=14;function kh(){}
_=kh.prototype=new Em();_.tI=0;_.a=null;function ch(b,a){b.a=a;return b}
function bh(){}
_=bh.prototype=new kh();_.tI=0;function ph(b,a){b.b=a;return b}
function rh(a){if(!a.a){a.a=$doc.createElement(ib);ag(a.b.d,a.a,0);a.a.appendChild($doc.createElement(jb))}}
function oh(){}
_=oh.prototype=new Em();_.tI=0;_.a=null;_.b=null;function Bh(a){a.a=uq(new tq());return a}
function Dh(d,b){var c,a;c=(a=b[kb],a==null?-1:a);if(c<0){return null}return Fd(xq(d.a,c),5)}
function Eh(d,b){var c,a;c=(a=b[kb],a==null?-1:a);b[kb]=null;zq(d.a,c,null)}
function sh(){}
_=sh.prototype=new Em();_.tI=0;function vh(b,a){b.b=a;xh(b);return b}
function xh(a){while(++a.a<a.b.a.b){if(xq(a.b.a,a.a)!=null){return}}}
function yh(){return this.a<this.b.a.b}
function zh(){var a;if(this.a>=this.b.a.b){throw new zr()}a=Fd(xq(this.b.a,this.a),5);xh(this);return a}
function th(){}
_=th.prototype=new Em();_.u=yh;_.x=zh;_.tI=0;_.a=-1;_.b=null;function gj(){gj=as;kj=dr(new cr());lj=hr(new gr())}
function fj(b,a){gj();b.a=Dj(new wj());b.i=a;fk(b);return b}
function hj(){var b,a;gj();var c,d;for(d=(b=co(new bo(),mq(lj.a).b.a),Cp(new Bp(),b));np(d.a.a);){c=Fd((a=Fd(op(d.a.a),12),a.q()),5);if(c.v()){c.A()}}}
function jj(b){gj();var a,c;c=Fd(ap(kj,b),6);if(c){return c}a=null;if(b!=null){if(!(a=$doc.getElementById(b))){return null}}if(kj.d==0){lf(new aj())}if(!a){c=dj(new cj())}else{c=fj(new Fi(),a)}gp(kj,b,c);ir(lj,c);return c}
function Fi(){}
_=Fi.prototype=new ng();_.tI=15;var kj,lj;function aj(){}
_=aj.prototype=new Em();_.tI=16;function ej(){ej=as;gj()}
function dj(a){ej();fj(a,$doc.body);return a}
function cj(){}
_=cj.prototype=new Fi();_.tI=17;function Bn(a,b){var c;while(a.u()){c=a.x();if(b==null?c==null:mc(b,c)){return a}}return null}
function Dn(a){throw new xn()}
function En(b){var a;a=Bn(this.w(),b);return !!a}
function An(){}
_=An.prototype=new Em();_.k=Dn;_.l=En;_.tI=0;function tp(a){this.j(this.ab(),a);return true}
function sp(b,a){throw new xn()}
function up(a,b){if(a<0||a>=b){xp(a,b)}}
function vp(e){var a,b,c,d,f;if((e==null?null:e)===(this==null?null:this)){return true}if(!(e!=null&&Ed(e.tI,13))){return false}f=Fd(e,13);if(this.ab()!=f.ab()){return false}c=lp(new jp(),this);d=f.w();while(c.a<c.b.ab()){a=op(c);b=op(d);if(!(a==null?b==null:mc(a,b))){return false}}return true}
function wp(){var a,b,c;b=1;a=lp(new jp(),this);while(a.a<a.b.ab()){c=op(a);b=31*b+(c==null?0:oc(c));b=~~b}return b}
function xp(a,b){throw qm(new pm(),lb+a+mb+b)}
function yp(){return lp(new jp(),this)}
function ip(){}
_=ip.prototype=new An();_.k=tp;_.j=sp;_.eQ=vp;_.hC=wp;_.w=yp;_.tI=18;function uq(a){a.a=vd(qe,0,0,0,0);a.b=0;return a}
function vq(b,a){yd(b.a,b.b++,a);return true}
function xq(b,a){up(a,b.b);return b.a[a]}
function yq(c,b,a){for(;a<c.b;++a){if(Fr(b,c.a[a])){return a}}return -1}
function zq(d,a,b){var c;c=(up(a,d.b),d.a[a]);yd(d.a,a,b);return c}
function Bq(a){return yd(this.a,this.b++,a),true}
function Aq(a,b){if(a<0||a>this.b){xp(a,this.b)}this.a.splice(a,0,b);++this.b}
function Cq(a){return yq(this,a,0)!=-1}
function Dq(a){return up(a,this.b),this.a[a]}
function Eq(){return this.b}
function tq(){}
_=tq.prototype=new ip();_.k=Bq;_.j=Aq;_.l=Cq;_.t=Dq;_.ab=Eq;_.tI=19;_.a=null;_.b=0;function oj(a){a.a=vd(qe,0,0,0,0);a.b=0;return a}
function qj(e,d,a){var b,c;for(c=lp(new jp(),e);c.a<c.b.ab();){b=Fd(op(c),7);b.a.b=Am(a);b.a.c=Am(d-1);$wnd.google.visualization.events.trigger(b.a.d,ob,null)}}
function nj(){}
_=nj.prototype=new tq();_.tI=20;function Dj(a){a.a=vd(oe,0,5,4,0);return a}
function Ej(a,b){bk(a,b,a.b)}
function ak(b,c){var a;for(a=0;a<b.b;++a){if(b.a[a]==c){return a}}return -1}
function bk(d,e,a){var b,c;if(a<0||a>d.b){throw new pm()}if(d.b==d.a.length){c=vd(oe,0,5,d.a.length*2,0);for(b=0;b<d.a.length;++b){yd(c,b,d.a[b])}d.a=c}++d.b;for(b=d.b-1;b>a;--b){yd(d.a,b,d.a[b-1])}yd(d.a,a,e)}
function ck(c,b){var a;if(b<0||b>=c.b){throw new pm()}--c.b;for(a=b;a<c.b;++a){yd(c.a,a,c.a[a+1])}yd(c.a,c.b,null)}
function dk(b,c){var a;a=ak(b,c);if(a==-1){throw new zr()}ck(b,a)}
function wj(){}
_=wj.prototype=new Em();_.tI=0;_.a=null;_.b=0;function zj(b,a){b.b=a;return b}
function Bj(){return this.a<this.b.b-1}
function Cj(){if(this.a>=this.b.b){throw new zr()}return this.b.a[++this.a]}
function xj(){}
_=xj.prototype=new Em();_.u=Bj;_.x=Cj;_.tI=0;_.a=-1;_.b=null;function Ck(b,c,a){var d;d=wl(new il());d.d=c;og(jj(a.id),d);if(d){Ek(c)}return d}
function Ek(b){b.getSelection=function(){return this.gwt_vis.r()};b.setSelection=function(a){this.gwt_vis.E(a)}}
function Fk(d,c){$wnd[d]=function(a){this.gwt_vis=Ck(c,this,a)};$wnd[d].prototype.draw=function(a,b){this.gwt_vis.o(a,b)}}
function zk(){}
_=zk.prototype=new zg();_.tI=21;_.d=null;function dl(b){var a,c;c=[];for(a=0;a<b.length;++a){c[a]=b[a]}c.constructor=$wnd.Array;return c}
function wl(a){a.a=eh(new ah());Bg(a,a.a);return a}
function yl(b,c){var a,d;gi(this.a);if(c){this.a.i.style[pb]=c.backgroundColor||qb}for(a=0;a<b.getNumberOfColumns();++a){si(this.a,0,a,rb+b.getColumnLabel(a)+sb)}for(d=0;d<b.getNumberOfRows();++d){for(a=0;a<b.getNumberOfColumns();++a){ti(this.a,d+1,a,b.getFormattedValue(d,a))}}di(this.a,kl(new jl(),this))}
function zl(){return dl(wd(ne,0,-1,[{row:this.c.a,column:this.b.a}]))}
function Al(a){$wnd.alert(tb)}
function il(){}
_=il.prototype=new zk();_.o=yl;_.r=zl;_.E=Al;_.tI=22;_.b=null;_.c=null;function kl(b,a){b.a=a;return b}
function jl(){}
_=jl.prototype=new Em();_.tI=23;_.a=null;function ul(a){if($wnd.onLoadCallback!=undefined){$wnd.onLoadCallback(a)}}
function rl(){}
_=rl.prototype=new Em();_.tI=0;function vn(){}
_=vn.prototype=new Em();_.tI=3;function km(){}
_=km.prototype=new vn();_.tI=4;function cn(){}
_=cn.prototype=new km();_.tI=5;function Dl(){}
_=Dl.prototype=new cn();_.tI=25;function em(c,a){var b;b=new am();return b}
function am(){}
_=am.prototype=new Em();_.tI=0;function bm(){}
_=bm.prototype=new cn();_.tI=28;function mm(){}
_=mm.prototype=new cn();_.tI=30;function qm(b,a){return b}
function pm(){}
_=pm.prototype=new cn();_.tI=31;function Cm(){}
_=Cm.prototype=new Em();_.tI=29;function wm(a,b){a.a=b;return a}
function ym(a){return a!=null&&Ed(a.tI,11)&&Fd(a,11).a==this.a}
function zm(){return this.a}
function Am(a){var b,c;if(a>-129&&a<128){b=a+128;c=(um(),vm)[b];if(!c){c=vm[b]=wm(new sm(),a)}return c}return wm(new sm(),a)}
function sm(){}
_=sm.prototype=new Cm();_.eQ=ym;_.hC=zm;_.tI=32;_.a=0;function um(){um=as;vm=vd(pe,0,11,256,0)}
var vm;function qn(b,a){if(a==null)return false;return b==a||b.toLowerCase()==a.toLowerCase()}
function tn(a){if(!(a!=null&&Ed(a.tI,1))){return false}return String(this)==a}
function un(){return mn(this)}
_=String.prototype;_.eQ=tn;_.hC=un;_.tI=2;function gn(){gn=as;hn={};ln={}}
function jn(e){var a,b,c,d;d=e.length;c=d<64?1:~~(d/32);a=0;for(b=0;b<d;b+=c){a<<=1;a+=e.charCodeAt(b)}a|=0;return a}
function mn(c){gn();var a=ub+c;var b=ln[a];if(b!=null){return b}b=hn[a];if(b==null){b=jn(c)}nn();return ln[a]=b}
function nn(){if(kn==256){hn=ln;ln={};kn=0}++kn}
var hn,kn=0,ln;function xn(){}
_=xn.prototype=new cn();_.tI=35;function mq(b){var a;a=ho(new ao(),b);return bq(new Ap(),b,a)}
function nq(c){var a,b,d,e,f;if((c==null?null:c)===(this==null?null:this)){return true}if(!(c!=null&&Ed(c.tI,14))){return false}e=Fd(c,14);if(Fd(this,14).d!=e.d){return false}for(b=co(new bo(),ho(new ao(),e).a);np(b.a);){a=Fd(op(b.a),12);d=a.q();f=a.s();if(!(d==null?Fd(this,14).c:d!=null&&Ed(d.tI,1)?cp(Fd(this,14),Fd(d,1)):bp(Fd(this,14),d,~~oc(d)))){return false}if(!Fr(f,d==null?Fd(this,14).b:d!=null&&Ed(d.tI,1)?Fd(this,14).e[ub+Fd(d,1)]:Eo(Fd(this,14),d,~~oc(d)))){return false}}return true}
function oq(){var a,b,c;c=0;for(b=co(new bo(),ho(new ao(),Fd(this,14)).a);np(b.a);){a=Fd(op(b.a),12);c+=a.hC();c=~~c}return c}
function zp(){}
_=zp.prototype=new Em();_.eQ=nq;_.hC=oq;_.tI=0;function zo(g,c){var e=g.a;for(var d in e){if(d==parseInt(d)){var a=e[d];for(var f=0,b=a.length;f<b;++f){c.k(a[f])}}}}
function Ao(e,a){var d=e.e;for(var c in d){if(c.charCodeAt(0)==58){var b=xo(e,c.substring(1));a.k(b)}}}
function Do(b,a){return a==null?b.c:a!=null&&Ed(a.tI,1)?cp(b,Fd(a,1)):bp(b,a,~~oc(a))}
function ap(b,a){return a==null?b.b:a!=null&&Ed(a.tI,1)?b.e[ub+Fd(a,1)]:Eo(b,a,~~oc(a))}
function Eo(h,g,e){var a=h.a[e];if(a){for(var f=0,b=a.length;f<b;++f){var c=a[f];var d=c.q();if(h.p(g,d)){return c.s()}}}return null}
function bp(h,g,e){var a=h.a[e];if(a){for(var f=0,b=a.length;f<b;++f){var c=a[f];var d=c.q();if(h.p(g,d)){return true}}}return false}
function cp(b,a){return ub+a in b.e}
function gp(b,a,c){return a==null?ep(b,c):a!=null&&Ed(a.tI,1)?fp(b,Fd(a,1),c):dp(b,a,c,~~oc(a))}
function dp(i,g,j,e){var a=i.a[e];if(a){for(var f=0,b=a.length;f<b;++f){var c=a[f];var d=c.q();if(i.p(g,d)){var h=c.s();c.F(j);return h}}}else{a=i.a[e]=[]}var c=sr(new rr(),g,j);a.push(c);++i.d;return null}
function ep(b,c){var a;a=b.b;b.b=c;if(!b.c){b.c=true;++b.d}return a}
function fp(d,a,e){var b,c=d.e;a=ub+a;if(a in c){b=c[a]}else{++d.d}c[a]=e;return b}
function hp(a,b){return (a==null?null:a)===(b==null?null:b)||a!=null&&mc(a,b)}
function Fn(){}
_=Fn.prototype=new zp();_.p=hp;_.tI=0;_.a=null;_.b=null;_.c=false;_.d=0;_.e=null;function rq(b){var a,c,d;if((b==null?null:b)===(this==null?null:this)){return true}if(!(b!=null&&Ed(b.tI,15))){return false}c=Fd(b,15);if(c.ab()!=this.ab()){return false}for(a=c.w();a.u();){d=a.x();if(!this.l(d)){return false}}return true}
function sq(){var a,b,c;a=0;for(b=this.w();b.u();){c=b.x();if(c!=null){a+=oc(c);a=~~a}}return a}
function pq(){}
_=pq.prototype=new An();_.eQ=rq;_.hC=sq;_.tI=36;function ho(b,a){b.a=a;return b}
function jo(c){var a,b,d;if(c!=null&&Ed(c.tI,12)){a=Fd(c,12);b=a.q();if(Do(this.a,b)){d=ap(this.a,b);return fr(a.s(),d)}}return false}
function ko(){return co(new bo(),this.a)}
function lo(){return this.a.d}
function ao(){}
_=ao.prototype=new pq();_.l=jo;_.w=ko;_.ab=lo;_.tI=37;_.a=null;function co(c,b){var a;c.b=b;a=uq(new tq());if(c.b.c){vq(a,no(new mo(),c.b))}Ao(c.b,a);zo(c.b,a);c.a=lp(new jp(),a);return c}
function fo(){return np(this.a)}
function go(){return Fd(op(this.a),12)}
function bo(){}
_=bo.prototype=new Em();_.u=fo;_.x=go;_.tI=0;_.a=null;_.b=null;function jq(b){var a;if(b!=null&&Ed(b.tI,12)){a=Fd(b,12);if(Fr(this.q(),a.q())&&Fr(this.s(),a.s())){return true}}return false}
function kq(){var a,b;a=0;b=0;if(this.q()!=null){a=oc(this.q())}if(this.s()!=null){b=oc(this.s())}return a^b}
function hq(){}
_=hq.prototype=new Em();_.eQ=jq;_.hC=kq;_.tI=38;function no(b,a){b.a=a;return b}
function po(){return null}
function qo(){return this.a.b}
function ro(a){return ep(this.a,a)}
function mo(){}
_=mo.prototype=new hq();_.q=po;_.s=qo;_.F=ro;_.tI=39;_.a=null;function to(c,a,b){c.b=b;c.a=a;return c}
function vo(){return this.a}
function wo(){return this.b.e[ub+this.a]}
function xo(b,a){return to(new so(),a,b)}
function yo(a){return fp(this.b,this.a,a)}
function so(){}
_=so.prototype=new hq();_.q=vo;_.s=wo;_.F=yo;_.tI=40;_.a=null;_.b=null;function lp(b,a){b.b=a;return b}
function np(a){return a.a<a.b.ab()}
function op(a){if(a.a>=a.b.ab()){throw new zr()}return a.b.t(a.a++)}
function pp(){return this.a<this.b.ab()}
function qp(){return op(this)}
function jp(){}
_=jp.prototype=new Em();_.u=pp;_.x=qp;_.tI=0;_.a=0;_.b=null;function bq(b,a,c){b.a=a;b.b=c;return b}
function eq(a){return Do(this.a,a)}
function fq(){var a;return a=co(new bo(),this.b.a),Cp(new Bp(),a)}
function gq(){return this.b.a.d}
function Ap(){}
_=Ap.prototype=new pq();_.l=eq;_.w=fq;_.ab=gq;_.tI=41;_.a=null;_.b=null;function Cp(a,b){a.a=b;return a}
function Fp(){return np(this.a.a)}
function aq(){var a;return a=Fd(op(this.a.a),12),a.q()}
function Bp(){}
_=Bp.prototype=new Em();_.u=Fp;_.x=aq;_.tI=0;_.a=null;function dr(a){a.a=[];a.e={};a.c=false;a.b=null;a.d=0;return a}
function fr(a,b){return (a==null?null:a)===(b==null?null:b)||a!=null&&mc(a,b)}
function cr(){}
_=cr.prototype=new Fn();_.tI=42;function hr(a){a.a=dr(new cr());return a}
function ir(c,a){var b;b=gp(c.a,a,c);return b==null}
function kr(b){var a;return a=gp(this.a,b,this),a==null}
function lr(a){return Do(this.a,a)}
function mr(){var a;return a=co(new bo(),mq(this.a).b.a),Cp(new Bp(),a)}
function nr(){return this.a.d}
function gr(){}
_=gr.prototype=new pq();_.k=kr;_.l=lr;_.w=mr;_.ab=nr;_.tI=43;_.a=null;function sr(b,a,c){b.a=a;b.b=c;return b}
function ur(){return this.a}
function vr(){return this.b}
function xr(b){var a;a=this.b;this.b=b;return a}
function rr(){}
_=rr.prototype=new hq();_.q=ur;_.s=vr;_.F=xr;_.tI=44;_.a=null;_.b=null;function zr(){}
_=zr.prototype=new cn();_.tI=45;function Fr(a,b){return (a==null?null:a)===(b==null?null:b)||a!=null&&mc(a,b)}
function Bl(){!!$stats&&$stats({moduleName:$moduleName,subSystem:vb,evtGroup:wb,millis:(new Date()).getTime(),type:xb,className:zb});Fk(Ab,new rl());ul($moduleName)}
function gwtOnLoad(b,d,c){$moduleName=d;$moduleBase=c;if(b)try{Bl()}catch(a){b(d)}else{Bl()}}
function as(){}
var oe=em(Bb,Cb),pe=em(Db,Eb),qe=em(Db,Fb),ne=em(ac,bc);$stats && $stats({moduleName:'com.google.gwt.visualization.customvisualization.CustomVisualization',subSystem:'startup',evtGroup:'moduleStartup',millis:(new Date()).getTime(),type:'moduleEvalEnd'});if (com_google_gwt_visualization_customvisualization_CustomVisualization) {  var __gwt_initHandlers = com_google_gwt_visualization_customvisualization_CustomVisualization.__gwt_initHandlers;  com_google_gwt_visualization_customvisualization_CustomVisualization.onScriptLoad(gwtOnLoad);}})();