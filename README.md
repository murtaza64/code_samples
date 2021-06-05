# Murtaza's Code Samples / portfolio
On this page you will find a subset of my programming skills, and hopefully what you're looking for is here.

The samples from Essaystance are less interesting than the snippets from random projects and challenges I've collected over the years, which often show more cleverness (perhaps at the cost of readability).

&nbsp;

[Coding Challenges](#coding-challenges)
Some selected coding challenge solutions I had lying around. Show problem solving, recursion, and use of Python-specific features.

[Essaystance (2020)](#essaystance)
My latest web project in Django. Snippets show some backend view logic and frontend.

[chess.py (2017)](#chess)
Terminal chess I built in stock Python. Uses OOP and generators extensively.

[RPi Timetable (2016)](#rpi-timetable)
Snippet shows a custom data structure used to reprsent an LED grid for displaying text. Employs multithreading.

[Friendship Graphs (2018)](#friendship-graphs)
A personal dataviz experiment/project. Snippet shows data manipulation and cleaning.

&nbsp;

Academic:

[Genetic Algorithm (Python, 2021)](#genetic-algorithm)
For CS470 AI, A genetic algorithm built from scratch for the Santa Fe/John Muir artificial ant problem

[Ray Tracing (C++, 2020)](#ray-tracing)
For CS478 Graphics, a simple raytracer in stock C++.

[Fork (C, 2020)](#fork)
For CS323 Systems, kernel code to handle forking a process.

[Hash Table (C, 2018)](#hash-table)
For CS223. Shows use of pointers and understanding of data structures.

&nbsp;

External links:

[Teleport (Python, 2018)](https://github.com/murtaza64/teleport/blob/master/tp.py)
A command line utility that allows you to teleport around your filetree! (full code linked)

The full source for [codeline](http://thecodeline.co), an older project of mine, is [up here](http://github.com/murtaza64/codeline). 

Some selected source files are included in this repo.

&nbsp;

Notes: 
* I forgot to include `git` on some versions of my resume. I indeed have experience with `git`, and have used it extensively in my projects. 
* The gap during 2019 was for medical reasons.
* I have experience with natural language processing, but so far have not done a substantial project requiring me to write more than small amounts of code.


## Coding Challenges
These are some snippets that I was pretty proud of when I wrote them.
### Nested Brackets (2017)
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

### Matrix determinant (2017)

```python
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

### Infix to postfix (2017)
Uses recursion
```python
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

### Piglatin oneliner (2016)
Legit use of a generator expression!
```python
def pig_it(text):
    return ' '.join(w[1:]+w[0]+'ay' if w.isalpha() else w for w in text.split(' '))
```

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
    doc = gapi_service.documents().get(documentId=doc_id).execute() 
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
 Handles the frontend rendering of the 1-5 star rater, allowing you to select a rating. When you hover over or click a star, all lower stars get filled in (so it looks pretty!). (The stars get darker when you click one to select a rating.) This snippet demonstrates the ability to use DOM event handlers to create a nice custom UI element.
 ```javascript
 function change_rating(id, value) { //this handler is attached in the HTML
    $(".star_" + id).each(function(index, item) {
        if (index < value) {
            $(item).attr("src", star_gold_url); //we change the star color by changing the source image
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
I created this terminal chess a while back. It makes use of some nice Python features. It supports everything except en passant and castling, which slipped my mind at the time. Below are some snippets, but the full code is available in this repo.

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

A few years ago, I built an LED timetable for my school desk using a Raspberry Pi. Rendering of text to the LED grid was handled by the classes in `matrix.py`, and the main script ran two threads (one for the modulation of the LEDs and one for updating the current lesson based on the current time). The [full script](https://github.com/murtaza64/timetable) includes a helper function to make a bitmap which I used for testing (and fun). It's mostly matrix transposition, but it was interesting coming up with a data structure and relevant methods to represent a scrolling LED sign. This is a little older, so please excuse the inadherence to Python style.

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
	
#...

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

	#... (code omitted to keep it short)
	
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

```

## Genetic Algorithm

This was a genetic algorithm I built for the artificial ant problem, wherein an ant with one perceptor and a few states is evolved to follow trails of food effectively. I have shown three functions here which are responsible for generating the ants of the next generation (the full `geneticAlgorithm.py` file is in the repo). `copulate` produces offspring between two parents based on a crossover probability. `mutate` changes genes of one ant based on a mutation probability. `select_parents` chooses parents from the set of ants in the current generation via a fitness rank lottery.

```python
def copulate(probability_crossover, *parents):
    '''
    generate a child from multiple parents
    '''
    parent_index = random.randrange(len(parents))
    child = ''
    for i in range(30):
        child += parents[parent_index][i]
        if random.random() < probability_crossover:
            parent_index = (parent_index + 1) % len(parents)
    return child

def mutate(ants, probability_mutation):
    '''
    ants : list[str] -- ants to mutate
    probability_mutation : number -- how likely each character is to randomly mutate

    '''
    new_ants = []
    mutations = 0
    for ant in ants:
        new_ant = ''
        for i, c in enumerate(ant):
            if random.random() < probability_mutation:
                if i % 3 == 0:
                    new_ant += str(random.randint(1, 4))
                else:
                    new_ant += str(random.randint(0, 9))
            else:
                new_ant += c
        if ant != new_ant:
            mutations += 1
        new_ants.append(new_ant)
    return new_ants, mutations

def select_parents(old_gen, n_parents, auto_top_n=-1):
    '''
    old_gen : list[str] -- population to select from
    n_parents : int -- how many parents to select

    return : list[str] -- population of ants selected
    '''
    if auto_top_n == -1:
        auto_top_n = n_parents//4

    #preselect the top n performers (default is a quarter of the parents to be selected)
    chosen_indices = set(range(auto_top_n))
    indices = list(range(len(old_gen)))

    #lottery for the rest
    #higher weights for higher fitness
    weights = [(population_size - i)**2  for i in indices]
    for index in range(auto_top_n):
        weights[index] = 0

    remaining = n_parents
    while len(chosen_indices) != n_parents:
        remaining = n_parents - len(chosen_indices)
        for selected_index in random.choices(indices, weights, k=remaining):
            if weights[selected_index]:
                weights[selected_index] = 0
                chosen_indices.add(selected_index)


    return [old_gen[i] for i in chosen_indices]
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
## Fork 
As part of CS323, we had to implement a bunch of kernel functionality on a toy OS. This is how I implemented the `fork` system call.
```C
case INT_SYS_FORK: {
        int found_slot = 0;
        int i;
        for (i = 1; i < NPROC; i++){
            if (processes[i].p_state == P_FREE) {
                found_slot = 1;
                break;
            }
        }
        if (!found_slot){
            current->p_registers.reg_rax = -1;
        }
        pid_t new_pid = processes[i].p_pid;
        proc* new = &processes[i];
        int failed = 0;

        //begin copied code [NOTE: this comment was for me, and was referring to my own code that I copied from earlier in the program]
        process_init(&processes[new_pid], 0);
        uintptr_t page_table_addrs[5];
        x86_64_pagetable* proc_pagetable_ptrs[5];
        for(i = 0; i < 5; i++){
            page_table_addrs[i] = get_available_page();
            if (page_table_addrs[i] == 0) {
                failed = 1;
                break;
            }
            assign_physical_page((uintptr_t) page_table_addrs[i], new_pid);
            proc_pagetable_ptrs[i] = (x86_64_pagetable*) page_table_addrs[i];
            memset(proc_pagetable_ptrs[i], 0, PAGESIZE);
        }
        
        x86_64_pagetable* proc_pagetable = proc_pagetable_ptrs[0];
        proc_pagetable_ptrs[0]->entry[0] =
            (x86_64_pageentry_t) proc_pagetable_ptrs[1] | PTE_P | PTE_W | PTE_U;
        proc_pagetable_ptrs[1]->entry[0] =
            (x86_64_pageentry_t) proc_pagetable_ptrs[2] | PTE_P | PTE_W | PTE_U;
        proc_pagetable_ptrs[2]->entry[0] =
            (x86_64_pageentry_t) proc_pagetable_ptrs[3] | PTE_P | PTE_W | PTE_U;
        proc_pagetable_ptrs[2]->entry[1] =
            (x86_64_pageentry_t) proc_pagetable_ptrs[4] | PTE_P | PTE_W | PTE_U;

        log_printf("address of pt for pid %i is %p\n", new_pid, proc_pagetable);
        log_printf("owner of pt is %i\n", pageinfo[PAGENUMBER(proc_pagetable)].owner);

        for(int addr = 0x0; addr < PROC_START_ADDR; addr += PAGESIZE) {
            vamapping map = virtual_memory_lookup(kernel_pagetable, addr);
            virtual_memory_map(proc_pagetable, addr, map.pa, PAGESIZE, map.perm, NULL);
        }

        processes[new_pid].p_pagetable = proc_pagetable;
        
        //++pageinfo[PAGENUMBER(proc_pagetable)].refcount;
        // int r = program_load(&processes[pid], program_number, NULL);
        // assert(r >= 0);
        processes[new_pid].p_registers.reg_rsp = MEMSIZE_VIRTUAL;
        uintptr_t stack_pa = get_available_page();
        uintptr_t stack_page = processes[new_pid].p_registers.reg_rsp - PAGESIZE;
        assign_physical_page(stack_pa, new_pid);
        virtual_memory_map(processes[new_pid].p_pagetable, stack_page, stack_pa,
                        PAGESIZE, PTE_P | PTE_W | PTE_U, NULL);
        processes[new_pid].p_state = P_RUNNABLE;
        //end copied code
	
        for (uintptr_t va = PROC_START_ADDR; va < MEMSIZE_VIRTUAL; va += PAGESIZE) {
            vamapping map = virtual_memory_lookup(current->p_pagetable, va);
            if (map.pn >= 0) {
                if ((map.perm & PTE_W) == 0) {
                    virtual_memory_map(processes[new_pid].p_pagetable, va, map.pa, PAGESIZE, map.perm, NULL);
                    pageinfo[PAGENUMBER(map.pa)].refcount++;
                }
                else {
                    uintptr_t new_pa = get_available_page();
                    if (new_pa == 0) {
                        failed = 1;
                        break;
                    }
                    assign_physical_page(new_pa, new_pid);
                    virtual_memory_map(processes[new_pid].p_pagetable, va, new_pa, PAGESIZE, map.perm, NULL);
                    memcpy((void*) new_pa, (void*) map.pa, PAGESIZE);
                }
            }
        }
        if (failed) {
            log_printf("failed a fork of pid %i\n", current->p_pid);
            current->p_registers.reg_rax = -1;
            for (uintptr_t pa = PAGESIZE; pa < MEMSIZE_PHYSICAL; pa += PAGESIZE) {
                // vamapping map = virtual_memory_lookup(new->p_pagetable, va);
                int pn = PAGENUMBER(pa);
                if (pageinfo[pn].owner == new->p_pid) {
                    log_printf("freeing page %i\n", pn);
                    pageinfo[pn].owner = 0;
                    pageinfo[pn].refcount = 0;
                }
            }
            new->p_state = P_FREE;
            break;
        }
        current->p_registers.reg_rax = new_pid;
        new->p_registers = current->p_registers;
        new->p_registers.reg_rax = 0;
        break;
```

## Friendship Graphs
I set up a little experiment to try and visualize the friend groups in our school. [This](https://i.imgur.com/wQLav5o.png) is how it turned out. I sent out a survey and used Gephi to visualize the resulting data. The Python to process the data using `networkx` was as follows:

```python
import json
import csv
from collections import namedtuple
import networkx as nx

G = nx.Graph()

with open('fixednames.json') as f:
    jnames = json.loads(f.read())
Name = namedtuple('Name', 'first last form year gender realFirst realLast uname')
Name.__str__ = lambda self: self.uname
names = [Name(**name) for name in jnames if name['year']=='12']
print(names)
#namesohjun = sorted([name for name in jnames if name['year']=='12'], key=lambda n: n['first'])

G.add_nodes_from(names)

with open('data.csv', newline='') as f:
    fr = csv.reader(f)
    header = next(fr)
    rows = reversed(list(fr))
    done_names = []
    for row in rows:
        date = row[0]
        thisname = row[1]
        if not thisname or thisname in done_names:
            continue
        else:
            done_names.append(thisname)
        for node in G.nodes():
            if node.first + ' ' + node.last == thisname.upper():
                thisnode = node
                print(thisnode)
                break
        else:
            raise ValueError('you messed up')
        for i, score in enumerate(row[2:]):
            if score == '':
                score = '0'
            try:
                G[thisnode][names[i]]['weight'] = (int(score) + G[thisnode][names[i]]['weight'])/2
            except KeyError:
                G.add_edge(thisnode, names[i], weight=int(score))

if __name__ == '__main__':
    # for node in nodes:
    #     node.print_friendships()

    for edge in G.edges(): #remove zeroes
        if G[edge[0]][edge[1]]['weight'] <= 2:
            G.remove_edge(*edge)
    for edge in G.edges(): #remove selflinks
        if edge[0] == edge[1]:
            G.remove_edge(*edge)

    for n, nbrs in G.adjacency_iter():
        #print(n.first, end=':\n')
        for nbr, eattr in nbrs.items():
            pass
            #print('\t', nbr.first, eattr['weight'])
    print(len(names))
    for name in names:
        for name2 in names:
            if str(name) == str(name2) and name != name2:
                print(name)
    nx.write_graphml(G, 'friendships.graphml')
```
## Hash Table
Here's part of an implementation of a hash table in C for CS 223 Data Structures and Programming Techniques. Nothing too special except demonstrating comfort with pointers. The full `smap.c` is included in the repo.
```c

bool smap_reallocate(smap* m) {
    int new_cap = m->arr_cap * 2;
    node** new_arr = smap_create_array(new_cap);
    if (new_arr == NULL) {
        return 0;
    }
    for (int i = 0; i < m->arr_cap; i++) {
        node* old_curr = m->arr[i];
        while(old_curr->next != NULL){
            int index = m->hash(old_curr->key) % new_cap;
            // printf("(%s=>%i) moved to m->arr[%i]\n", old_curr->key, *(old_curr->val), index);
            node* new_curr = new_arr[index];
            while (new_curr->next != NULL) {
                new_curr = new_curr->next;
            }
            new_curr->key = old_curr->key;
            new_curr->val = old_curr->val;
            new_curr->next = smap_create_empty_node();
            if (new_curr->next == NULL) {
                smap_destroy_array(new_arr, new_cap, 0);
                return 0;
            }
            old_curr = old_curr->next;
        }
    }
    smap_destroy_array(m->arr, m->arr_cap, 0);
    m->arr = new_arr;
    m->arr_cap = new_cap;
    return 1;
}

bool smap_put(smap *m, const char *key, int *value){
    // printf("PUTTING (%s=>%i)\n  ", key, *value);
    if (key == NULL || m == NULL) {
        return 0;
    }
    node* curr = smap_find_key_node(m, key);
    if (curr->key == NULL) {
        if (m->size+1 == m->arr_cap) {
            //reallocate array and rehash elements
            smap_reallocate(m);
            curr = smap_find_key_node(m, key);
        }
        curr->next = smap_create_empty_node();
        if (curr->next == NULL) {
            return 0;
        }
        char* new_key = malloc(sizeof(char) * (strlen(key)+1));
        if (new_key == NULL) {
            free(curr->next);
            curr->next = NULL;
            return 0;
        }
        strcpy(new_key, key);
        curr->key = new_key;
        m->size++;
    }
    curr->val = value;
    return 1;
}

int *smap_get(smap *m, const char *key) {
    if (key == NULL || m == NULL) {
        return NULL;
    }
    node* curr = smap_find_key_node(m, key);
    if (curr->key == NULL) {
        return NULL;
    }
    else {
        return curr->val;
    }
}
```

