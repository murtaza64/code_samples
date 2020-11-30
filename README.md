# Murtaza's Code Samples / portfolio
The samples from Essaystance are less interesting than the snippets from random projects and challenges I've collected over the years, which sometimes show more cleverness (at the cost of readability).


* [Essaystance](#essaystance)
* [chess.py](#chess)
* [RPi Timetable](#rpi-timetable)
* [Ray Tracing (C++)](#ray-tracing)
* [Coding Challenges](#coding-challenges)

Small note: I forgot to include `git` on my resume. I indeed have experience with `git`, and have used it extensively in my projects.

## Essaystance
This platform provides highschoolers with high quality, affordable essay feedback from students at the universities they want to apply to.

### Priority filter
Method of `DispatchView` to return a list of readers that match the criteria. Assigns each reader a score based on how relevant they are to the search, so that AND matches are better than OR matches, but both are listed.
```python
def dispatch_get_readers(request, pk):
    if request.method == 'POST':
        l = []
        fil = request.POST['filter']
        if not fil:
            readers = Reader.objects.all()
        else:
            readers = []
            reader_qs = Reader.objects.all()
            score = {}
            colleges = []
            majors = []
            fil_obj = json.loads(fil)
            
            #extract filter parameters from the request
            for param in fil_obj:
                if param['type'] == 'college':
                    try:
                        colleges.append(College.objects.get(id=param['id']))
                    except College.DoesNotExist:
                        pass
                elif param['type'] == 'major':
                    try:
                        majors.append(Major.objects.get(id=param['id']))
                    except Major.DoesNotExist:
                        pass
                        
            #assign each reader a score
            for reader in reader_qs:
                score[reader.id] = 0
                for college in colleges:
                    if reader.college == college:
                        score[reader.id] += 10
                for major in majors:
                    if major in reader.majors.all():
                        score[reader.id] += 4
                if score[reader.id] > 0: #only return readers that match at least one parameter
                    readers.append(reader)

            readers.sort(key=lambda r: score[r.id], reverse=True)

        essay = Essay.objects.get(id=pk)
        
        #only consider active readers
        for reader in filter(lambda r: r.active, readers):
            try:
                assert(EmailVerification.objects.get(user=reader.user).verified)
            except(AssertionError, EmailVerification.DoesNotExist):
                continue   
            obj = {
              first_name': reader.user.first_name
              last_name': reader.user.last_name
              rating': reader.rating
              'majors': [major.name for major in reader.majors.all()]
              graduating_year': reader.graduating_year
              college': reader.college.alt_name if reader.college.alt_name else reader.college.name
              id': reader.id
              price': f'{price_str(calculate_price(essay, reader))}'
            }
            l.append(obj)
        p = Paginator(l, 10)
        page_number = int(request.POST['page'])
        try:
            return JsonResponse({'success': True, 'objects': p.page(page_number).object_list})
        except EmptyPage:
            return JsonResponse({'success': False, 'code': 'empty_page'})
```
### gd_word_count_essay
Counts the number of words in a Google Doc, using Drive API and parsing the document manually (since the API provides no endpoint for word count)
```python
def gd_word_count_essay(essay, user):
    spl = essay.link.split('/')
    doc_id = spl[spl.index('d') + 1]
    credset = OAuthCredentialSet.objects.get(user=user) #OAuth credentials are stored in the Essaystance db
    credentials = credset.to_credentials()
    gapi_service = build('docs', 'v1', credentials=credentials)
    doc = gapi_service.documents().get(documentId=doc_id).execute() #TODO: error checking--timeout
    credset.update_from_credentials(credentials)
    n = 0
    for section in doc['body']['content']:
        try:
            for el in section['paragraph']['elements']:
                content = el['textRun']['content'].strip()
                if content:
                    n += len(content.split(' ')) #the meat of this function
        except KeyError:
            continue
    essay.wordcount = n
    gapi_service = build('drive', 'v3', credentials=credentials)
    permissions = gapi_service.permissions().list(fileId=doc_id).execute()
    viewable = False
    for perm in permissions['permissions']:
        if perm['type'] == 'anyone' and perm['role'] in ['reader', 'writer', 'owner']:
            viewable = True
            break
    if not viewable: #change the sharing rules if required
        gapi_service.permissions().create(fileId=doc_id, body={'role': 'reader', 'type': 'anyone'}).execute()

    essay.save()
 ```
 
 ### Rate your review (js)
 Handles the frontend rendering of the 1-5 star rater, allowing you to select a rating. When you hover over or click a star, all lower stars get filled in (so it looks pretty!). (The stars get darker when you click one to select a rating.)
 ```javascript
 function change_rating(id, value) { //this handler is attached in the HTML
    $(".star_" + id).each(function(index, item) {
        if (index < value) {
            $(item).attr("src", star_gold_url);
        }
        else {
            $(item).attr("src", star_grey_url);
        }
    });
    ratings[id] = value;
    $("#rate_" + id).removeClass("disabled")
}

function star_hover(event) {
    id = event.target.className.substring(10)
    // console.log("HOVERING " + id)
    hover = true
    $(".star_" + id).each(function(index, item) {
        if (hover) {
            $(item).attr("src", star_light_url);
            if (item == event.target) {
                hover = false
            }
        }
        else {
            $(item).attr("src", star_grey_url);
        }
    });
}

function star_unhover(event) {
    id = parseInt($(event.target).parents(".reviewer").children(".req_id").html());
    value = ratings[id];
    $(".star_" + id).each(function(index, item) {
        $(item).attr("src", star_grey_url);
    });
    if (value !== undefined) {
        $(".star_" + id).each(function(index, item) {
            if (index < value) {
                $(item).attr("src", star_gold_url);
            }
            else {
                $(item).attr("src", star_grey_url);
            }
        });
    }
}

function enable_actions(id) { //we disable actions (link clicks) when submitting the rating, so we need to be able to reenable them if the request fails
    button = $("#rate_" + id);
    button.empty();
    button.html("Rate review");
    button.css({
        "padding": "",
        "background-color": "",
        "border": ""
    });
    button.addClass("fake_link");
    button.attr("onclick", "submit_rating(" + id + ")")
    $(".star_" + id).each(function(index, item) {
        $(item).parent().addClass("fake_link");
        $(item).parent().attr("onclick", "change_rating(" + id + ", " + (index + 1) + ")");
    });
    attach_hover_handlers();
}

function submit_rating(id) {
    $("#error_" + id).css("display", "none");
    button = $("#rate_" + id);
    button.css({
        "padding": "3px",
        "width": button.css("width"),
        "height": button.css("height"),
        "background-color": "#ffffff",
        "border": "2px solid #ffc107"
    });
    button.html("");
    button.removeClass("fake_link");
    button.append("<img src='" + spinner_url + "' width=15px height=auto style='margin:auto; display:block;'>");
    $(".star_" + id).each(function(index, item) {
        $(item).parent().removeClass("fake_link");
        $(item).parent().removeAttr("onclick");
    });
    button.parent().removeClass("fake_link");
    button.parent().removeAttr("onclick");
    $(".star_" + id).each(function(index, item) {
        $(item).off("mouseenter mouseleave")
    });
    $.ajax({
        type: 'POST',
        url: '/rate/' + id,
        data: {rating: ratings[id]},
        beforeSend: function(xhr, settings) {
            xhr.setRequestHeader("X-CSRFToken", csrftoken);
        },
        success: function(data, status, jqXHR) {
            // console.log(data);
            if (data.success) {
                button.empty();
                button.css("text-align", "center");
                button.html("Thanks!")
            }
            else {
                console.log("server error rating req " + id);
                console.log(data.message);
                $("#error_" + id).css("display", "block");
                enable_actions(id);
            }
        },
        error: function(jqXHR, status, errorThrown) {
            console.log("failed to POST rating req " + id);
            console.log(status, errorThrown);
            $("#error_" + id).css("display", "block");
            enable_actions(id);
        }
    })
}

function attach_hover_handlers(){
    $(".reviewer.dashboard").each(function(index, item) {
        id = parseInt($(item).children(".req_id").html()) //we need to attach separate hover handlers for each review rater
        // console.log("reviewer " + id)
        $(item).children().children(".stars_container").on("mouseleave", star_unhover)
        $(".star_" + id).each(function(index2, item2) {
            // console.log(id);
            $(item2).on({
                mouseenter: star_hover,
            })
        });
    });
}

$(document).ready(function() {
    radio_click();
    attach_hover_handlers();
})
```

## Chess 
I created this terminal chess a while back. It makes use of some nice Python features. It supports everything except en passant and castling, which slipped my mind at the time. Below are some snippets, but the full code is available in the repo.

Example piece class. `theoretical moves` is a generator used to retrieve all moves that the piece could make if the board were empty. This is combined with `validate_move` to create the square highlights in the terminal when you select a piece. `validate_move` is sometimes used directly if algebraic chess notation is used to make the move.
```python
class Pawn(Piece):

    symbol = 'P'
    name = 'Pawn'

    def validate_move(self, x2, y2, board):
        try:
            x1, y1 = self.get_coords(board)
            assert(Piece.validate_move(self, x2, y2, board))
            vert_op = op.add if self.color == Side.BLACK else op.sub
            if board[x2,y2]:
                assert(x2 in (x1 - 1, x1 + 1))
                assert(y2 == vert_op(y1, 1))
            else:
                assert(x2 == x1)
                if y1 == (board.dim - 2 if self.color == Side.WHITE else 1):
                    assert(y2 in (vert_op(y1, 1), vert_op(y1, 2)))
                else:
                    assert(y2 == vert_op(y1, 1))
            return True
        except AssertionError:
            return False
    
    def theoretical_moves(self, board):
        x1, y1 = self.get_coords(board)
        vert_op = op.add if self.color == Side.BLACK else op.sub
        if y1 == (board.dim - 2 if self.color == Side.WHITE else 1):
            yield (x1, vert_op(y1, 2))
        yield (x1, vert_op(y1, 1))
        yield (x1 + 1, vert_op(y1, 1))
        yield (x1 - 1, vert_op(y1, 1))

```

Here is the method of the `Board` class that prints the board, highlighting available moves of the current piece, using terminal control characters to change foreground and background color:

```python
def print_piece_info(self, piece):
    os.system('cls' if os.name == 'nt' else 'clear')
    x1, y1 = piece.get_coords(self)
    ret = '   '
    for x in 'abcdefgh':
        ret += ' ' + x + ' '
    ret += '\n'
    bg = self.bg_iter()
    for y, row in enumerate(self.grid):
        ret += str(self.dim-y) + '  '
        for x, cell in enumerate(row):
            if (x, y) == (x1, y1):
                next(bg)
                bgstr = BG.YELLOW
            elif (x, y) in piece.theoretical_moves(board):
                next(bg)
                if board[x1, y1].validate_move(x, y, board):
                    bgstr = BG.GREEN
                else:
                    bgstr = BG.RED
            else:
                bgstr = next(bg)
            if cell:
                ret += cell.color.fg
                ret += bgstr + ' ' + cell.symbol + ' '
            else:
                ret += bgstr + '   '
        ret += CLEAR_FMT + '  ' + str(self.dim-y) + '\n'
    ret += '   '
    for x in 'abcdefgh':
        ret += ' ' + x + ' '
    ret += CLEAR_FMT
    print(ret)
    for line in self._log[-12:]:
        print(line)
```
Check out the rest of the script for more.

## RPi Timetable

A few years ago, I built an LED timetable for my school desk using a Raspberry Pi. Rendering of text to the LED grid was handled by the classes in `matrix.py`, and the main script ran two threads (one for the modulation of the LEDs and one for updating the current lesson based on the current time). The script includes a helper function to make a bitmap which I used for testing (and fun). This is a little older, so please excuse the inadherence to Python style.

```python
from time import sleep
from struct import *
from scrollalphabet import *
import sys, os

errors=''
EMPTYTRIPLE=(0,0,0)
WHITE=(255,255,255)
def monochrome(inMatrix):
	matrix=[]
	for row in inMatrix:
		matrix.append([])
		for cell in row:
			if cell==1:
				matrix[-1].append(WHITE)
			else:
				matrix[-1].append(EMPTYTRIPLE)
	return matrix
def makeBMP(matrix, name):
	global errors
	w, h = 0, 0
	for row in matrix:
		h+=1
	for cell in matrix[0]:
		w+=1
	depth = 24
	padding = int(4 - (w*(depth/8) %4)) if (w*(depth/8) % 4 != 0) else 0
	bytes_per_row = int (w*(depth/8) + padding) 
	bytes_in_image = h*bytes_per_row

	size = 54 + bytes_in_image
	BITMAPFILEHEADER = pack('<hihhi',0x4D42, size, 0, 0, 54)

	BITMAPINFOHEADER = pack('<iiihhiiiiii', 40, w, -h, 1, depth, 0, bytes_in_image, 2835, 2835, 0, 0)

	imagebytes = []

	for row in matrix:
		for cell in row:
			red=cell[0].to_bytes(1, byteorder='little')
			green=cell[2].to_bytes(1, byteorder='little')
			blue=cell[1].to_bytes(1, byteorder='little')

			imagebytes.append(blue)
			imagebytes.append(green)
			imagebytes.append(red)
			
		for i in range(0, padding):
			imagebytes.append((0).to_bytes(1, byteorder='little'))
	try:
		image=open('bmp/'+name, 'wb')
		image.write(BITMAPFILEHEADER)
		image.write(BITMAPINFOHEADER)
		for byte in imagebytes:
			image.write(byte)
		image.close()
	except PermissionError:
		print('permission denied in ' + 'bmp/' + name)
		errors += 'permission denied in ' + 'bmp/' + name + '\n'
	except FileNotFoundError:
		os.makedirs('bmp/')
	

class pixelmatrix:
	def __init__(self, dims=(8,8), BMP=0, _print=1):
		self.matrix=[]
		self.frame=0
		y=dims[0]
		x=dims[1]
		self.ylen=y
		self.xlen=x
		self.BMP=BMP
		self._print=_print
		for i in range (y):
			row=[EMPTYTRIPLE for k in range (x)]
			self.matrix.append(row)
	def __getitem__(self, key):
		return self.matrix[key]
	def step(self):
		if self._print:
			for row in self.matrix:
				print('[', end='')
				for cell in row:
					if cell == EMPTYTRIPLE:
						print(' ', end ='')
					else:
						print('X', end ='')
				print(']\n', end='')
			print('\n')
		filename='{0:0>3}.bmp'.format(self.frame)
		self.frame+=1
		sleep(0.1)
		if self.BMP:
			makeBMP(self.matrix, filename)

	def update(self, inMatrix, coords=(0,0)):
		y,x=coords
		ylen=len(inMatrix)
		xlen=len(inMatrix[0])
		for i, row in enumerate(inMatrix):
			for j, cell in enumerate(row):
				self.matrix[y+i][x+j] = cell

	def newRow(self):
		return [EMPTYTRIPLE for k in range (self.xlen)]
	def newCol(self):
		return [EMPTYTRIPLE for k in range (self.ylen)]
	def pushRow(self, row, side='bottom'):
		xin = len(row)
		if side == 'bottom':
			del self.matrix[0]
			self.matrix.append(self.newRow())
			for i in range (xin):
				self.matrix[-1][i]=row[i]
		if side == 'top':
			del self.matrix[-1]
			self.matrix.insert(0, self.newRow())
			for i in range (xin):
				self.matrix[0][i]=row[i]
	def pushCol(self, col, side='right'):
		yin = len(col)
		if side == 'right':
			for row in self.matrix:
				row.pop(0)
			for i in range (yin):
				self.matrix[i].append(col[i])
			for j in range (i+1, self.ylen):
				self.matrix[j].append(EMPTYTRIPLE)
		if side == 'left':
			for row in self.matrix:
				row.pop()
			for i in range (yin):
				self.matrix[i].insert(0, col[i])
			for j in range (i+1, self.ylen):
				self.matrix[j].insert(0, EMPTYTRIPLE)
	def scroll(self, inMatrix):
		width=len(inMatrix[0])
		for i in range (width):
			col=[]
			for row in inMatrix:
				col.append(row[i])
			self.pushCol(col)
			self.step()
			
	def scrollText(self, str):
		for char in str:
			chargrid=monochrome(scrollAlphabet[char])
			self.scroll(chargrid)
			self.pushCol(self.newCol())
			self.step()
		for i in range (self.xlen):
			self.pushCol(self.newCol())
			self.step()

if __name__ == "__main__":
	testStr='OO'

	#scrollMatrix(testStr)

	#updateMatrix(testMatrix, (4,4))
	#showMatrix()

	screen = pixelmatrix(dims=(8,16))
	screen.step()
	screen.scrollText(sys.argv[1])
	screen.update(monochrome(testMatrix), (4,4))
	screen.step()
	print(errors)
```

## Ray Tracing

This is from CS 478 Graphics. It generated [this movie](https://www.youtube.com/watch?v=UhzPQjyFjdE). The full file is available in the repo (`previz.cpp`), but this snippet shows how glass refraction is calculated for a sphere. This includes Fresnel effects.

```cpp
    else if (s.material == GLASS) {
      ray refract;
      refract.origin = p;
      VEC3 in = (p - r.origin)/(p - r.origin).norm();
      float costheta;
      float cosphi;
      float k_reflect;
      if (r.in_glass) {
        // printf("in glass\n");
        n = -n;
        // printf("%f\n", n[2]);
        VEC3 dir1 = (1.5/1) * (in - n * in.dot(n));
        VEC3 dir2 = sqrt(1 - pow(1.5/1, 2) * (1 - pow(in.dot(n), 2))) * n;
        refract.dir = (dir1 - dir2)/(dir1 - dir2).norm();
        refract.in_glass = 0;
        costheta = -in.dot(n);
        cosphi = sqrt(1 - pow((1.5/1 * sqrt(1 - costheta*costheta)), 2)); 
        k_reflect = fresnel(costheta, cosphi, 0);
      }
      else {
        // printf("hit glass\n");
        VEC3 dir1 = (1/1.5) * (in - n * in.dot(n));
        VEC3 dir2 = sqrt(1 - pow(1/1.5, 2) * (1 - pow(in.dot(n), 2))) * n;
        refract.dir = (dir1 - dir2)/(dir1 - dir2).norm();
        refract.in_glass = 1;
        costheta = -in.dot(n);
        cosphi = sqrt(1 - pow((1/1.5 * sqrt(1 - costheta*costheta)), 2));
        k_reflect = fresnel(costheta, cosphi, 1);
      }

      ray reflect;
      VEC3 e = (r.origin-p)/(r.origin-p).norm();
      reflect.origin = p;
      reflect.dir = (-e + 2*(n.dot(e))*n)/(-e + 2*(n.dot(e))*n).norm();


      VEC3 refract_color = rayColor(refract, nlights, lights, depth+1);
      VEC3 reflect_color = rayColor(reflect, nlights, lights, depth+1);
      color = k_reflect*reflect_color + (1-k_reflect)*refract_color;
    }
```

## Coding Challenges
### Nested Brackets
Confirm if a string's brackets are nested properly

```python
#st = '([{[]}][({}{})[]])'
b_map = {']':'[', '}':'{', ')':'(', '>':'<'}
def find_inners(current_bracket, l):
	inners = []
	while l:
		c = l.pop(0)
		#print(c)
		if c in '[{(':
			inners.append(find_inners(c, l))
			#print(inners)
		elif b_map[c] == current_bracket:
			return inners
		else: 
			raise SyntaxError
	raise SyntaxError

def check_brackets(st):
	try:
		find_inners('<', list(st)+['>'])
	except SyntaxError:
		return False
	return True
#print(check_brackets(st))
```

### Matrix determinant

```
def determinant(matrix):
    minor_dets = []
    for i, val in enumerate(matrix[0]):
        if len(matrix) <= 1:
            return matrix[0][0]
        minor_dets.append(val * determinant([row[:i]+row[i+1:] for row in matrix[1:]]))
    result = 0
    while minor_dets:
        try:
            result += minor_dets.pop(0)
            result -= minor_dets.pop(0)
        except IndexError:
            break
    return result
```

### Infix to postfix

```
def postfix_recursive(infix_l):
    ops = {'+': 1, '-': 1, '/': 2, '*': 2, '^': 3}
    postfix, cur_ops = [], []
    while infix_l:
        c = infix_l.pop(0)
        if c in ops:
            while cur_ops and ops[c] <= ops[cur_ops[-1]]:
                postfix.append(cur_ops.pop())
            cur_ops.append(c)
        elif c == '(':
            postfix += postfix_recursive(infix_l)
        elif c == ')':
            while cur_ops:
                postfix.append(cur_ops.pop())
            return postfix
        elif c.isdigit():
            postfix.append(c)

to_postfix = lambda infix: ''.join(postfix_recursive(list(infix+')')))
```

### Piglatin oneliner

```
def pig_it(text):
    return ' '.join(w[1:]+w[0]+'ay' if w.isalpha() else w for w in text.split(' '))
```
