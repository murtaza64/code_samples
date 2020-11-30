from enum import Enum 
from itertools import cycle, chain
import operator as op
import os
import random

class Side(Enum):
    WHITE = 1
    BLACK = 2

    @property
    def fg(self):
        if self == Side.WHITE:
            return '\x1b[37m'
        if self == Side.BLACK:
            return '\x1b[30m'

class BG:
    BLUE =   '\x1b[44m'
    CYAN =   '\x1b[46m'
    YELLOW = '\x1b[43m'
    GREEN =  '\x1b[42m'
    RED =    '\x1b[41m'

CLEAR_FMT = '\x1b[0m'

class ChessBoard():

    dim = 8

    class InvalidMove(Exception):
        pass

    def __init__(self, AI=False):
        self.grid = [[Rook(Side.BLACK), Knight(Side.BLACK), Bishop(Side.BLACK), Queen(Side.BLACK), King(Side.BLACK), Bishop(Side.BLACK), Knight(Side.BLACK), Rook(Side.BLACK)],
                     [Pawn(Side.BLACK) for x in range(self.dim)],
                     [None]*self.dim,
                     #[None, None, None, None, None, None, Queen(Side.WHITE), None],
                     #[None, None, None, None, King(Side.BLACK), None, None, None],
                     [None]*self.dim,
                     [None]*self.dim,
                     [None]*self.dim,
                     [Pawn(Side.WHITE) for x in range(self.dim)],
                     [Rook(Side.WHITE), Knight(Side.WHITE), Bishop(Side.WHITE), Queen(Side.WHITE), King(Side.WHITE), Bishop(Side.WHITE), Knight(Side.WHITE), Rook(Side.WHITE)]]
        self.turn_order = cycle(Side)
        self.turn = next(self.turn_order)
        self._log = []
        self.in_check = False
        self.AI = AI

    def __getitem__(self, indices):
        if isinstance(indices, int):
            return self.grid[indices//self.dim][indices%self.dim]
        elif isinstance(indices, tuple) and len(indices) == 2:
            try:
                return self.grid[indices[1]][indices[0]]
            except IndexError:
                return None
        else:
            raise TypeError('subscript must be a coordinate pair')

    def __setitem__(self, indices, item):
        if isinstance(indices, int):
            self.grid[indices//self.dim][indices%self.dim] = item
        elif isinstance(indices, tuple) and len(indices) == 2:
            self.grid[indices[1]][indices[0]] = item
        else:
            raise TypeError('subscript must be a coordinate pair or integer')
        
    @staticmethod
    def bg_iter():
        return cycle(chain([BG.CYAN, BG.BLUE]*4, [BG.BLUE, BG.CYAN]*4))

    def log(self, text):
        print(text)
        self._log.append(text)
    
    def get_input(self, prompt):
        print(prompt + ': ', end='')
        s = input()
        self._log.append(prompt + ': ' + s)
        return s
    
    def stringify_move(self, x1, y1, x2, y2):
        return 'abcdefgh'[x1] + str(self.dim - y1) + 'abcdefgh'[x2] + str(self.dim - y2)
    
    def stringify_coords(self, x, y):
        'abcdefgh'[x] + str(self.dim - y)

    def play(self):
        try:
            while True:
                self.print_board()
                try:
                    if self.AI:
                        if self.turn == Side.BLACK:
                            moves = [(*piece.get_coords(self), *move)
                                    for piece in self if piece and piece.color == Side.BLACK
                                    for move in piece.theoretical_moves(self) 
                                    if piece.validate_move(*move, self)]
                            mv = random.choice(moves)
                            self.log('AI (' + self.turn.name + '): ' + self.stringify_move(*mv))
                            self.move(*mv)
                            continue
                    coords = self.get_input('make a move (' + board.turn.name + ')')
                    message = self.algebraic(coords)
                    if isinstance(message, str):
                        self.log(message)
                    self.in_check = self.is_in_check()
                    if self.in_check:
                        if self.checkmate():
                            self.print_board()
                            self.log('checkmate!')
                            self.log((next(self.turn_order).name + ' wins!'))
                            break
                        else: 
                            self.log(self.turn.name + ' is in check!')
                except ChessBoard.InvalidMove as e:
                    self.log('invalid move: ' + str(e))
                except ValueError:
                    self.log('invalid move string')
        except (KeyboardInterrupt, EOFError):
            self.log('game ended')

    def is_in_check(self):
        for piece in self:
            if piece and isinstance(piece, King) and piece.color == self.turn:
                king = piece
        coords = king.get_coords(self)
        check = False
        for piece in self:
            if piece and piece.color != self.turn:
                for endcoords in piece.theoretical_moves(self):
                    if endcoords == coords and piece.validate_move(*endcoords, board):
                        check = True
        return check

    def checkmate(self):
        if not self.in_check:
            return False
        for piece in self:
            if piece and piece.color == self.turn:
                for endcoords in piece.theoretical_moves(self):
                    if piece.validate_move(*endcoords, board):
                        coords = piece.get_coords(self)
                        endpiece = self[endcoords]
                        self[endcoords] = piece
                        self[coords] = None
                        check = self.is_in_check()
                        self[endcoords] = endpiece
                        self[coords] = piece
                        if not check:
                            print(coords, endcoords)
                            return False
        return True


    def add_piece(self, piece, x, y):
        self.grid[y][x] = piece
    
    def algebraic(self, movestr):
        try:
            a = list(movestr)
            if len(a) == 3:
                sym = a[0]
                x2 = 'abcdefgh'.index(a[1])
                y2 = self.dim - int(a[2])
                move_choices = []
                for piece in self:
                    if piece and piece.color == self.turn and piece.symbol == sym.upper():
                        for move in piece.theoretical_moves(board):
                            if move == (x2, y2) and piece.validate_move(*move, board):
                                move_choices.append((piece, move))
                if not move_choices:
                    return 'no matching valid move'
                elif len(move_choices) == 1:
                    move = (*move_choices[0][0].get_coords(self), *move_choices[0][1])
                    self.move(*move)
                    return 'successful move ' + self.stringify_move(*move)
                else:
                    piece_locations = [i[0].get_coords(self) for i in move_choices]
                    piece_locations_strs = [self.stringify_coords(x, y) for x, y in piece_locations]
                    self.print_ambiguous_highlights(piece_locations, (x2, y2))
                    s = self.get_input('ambiguous move; choose from ' + sym + ' at ' + ', '.join(piece_locations_strs))
                    if s in piece_locations_strs:
                        self.move('abcdefgh'.index(s[0]), self.dim - int(s[1]), x2, y2)
                        return 'successful move ' + self.stringify_move(*move)
            elif len(a) == 4:
                x1 = 'abcdefgh'.index(a[0])
                x2 = 'abcdefgh'.index(a[2])
                y1 = self.dim - int(a[1])
                y2 = self.dim - int(a[3])
                self.move(x1, y1, x2, y2)
                return 'successful move ' + movestr
            elif len(a) == 2:
                x1 = 'abcdefgh'.index(a[0])
                y1 = self.dim - int(a[1])
                piece = self[x1, y1]
                if piece:
                    self.print_piece_info(piece)
                    if piece.color == self.turn:
                        movetostr = self.get_input('move this ' + str(self[x1, y1]) + ' to')
                        b = list(movetostr)
                        assert(len(b) == 2)
                        x2 = 'abcdefgh'.index(b[0])
                        y2 = self.dim - int(b[1])
                        self.move(x1, y1, x2, y2)
                        return 'successful move ' + movestr + movetostr
                    else:
                        self.get_input('wrong color. press enter to continue')
                else:
                    return 'empty space'
            else:
                raise ValueError
        except (ValueError, AssertionError):
            raise ValueError('invalid move string')

    def move(self, x1, y1, x2, y2):
        startpiece = self[x1, y1]
        endpiece = self[x2, y2]
        if not startpiece:
            raise ChessBoard.InvalidMove('no piece at start point')
        if not startpiece.color == self.turn:
            raise ChessBoard.InvalidMove('wrong color')
        if startpiece.validate_move(x2, y2, self):
            if endpiece:
                if endpiece.color != startpiece.color:
                    endpiece.dead = True
                else:
                    raise ChessBoard.InvalidMove('move not legal')
            self[x2, y2] = startpiece #change this to enable pawn promotion
            self[x1, y1] = None
            check = self.is_in_check()
            if check:
                self[x1, y1] = startpiece
                self[x2, y2] = endpiece
                raise ChessBoard.InvalidMove('move would endanger king')
            self.turn = next(self.turn_order)
        else:
            raise ChessBoard.InvalidMove('move not legal')

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
    
    def print_ambiguous_highlights(self, piece_locations, move):
        os.system('cls' if os.name == 'nt' else 'clear')
        ret = '   '
        for x in 'abcdefgh':
            ret += ' ' + x + ' '
        ret += '\n'
        bg = self.bg_iter()
        for y, row in enumerate(self.grid):
            ret += str(self.dim-y) + '  '
            for x, cell in enumerate(row):
                if (x, y) in piece_locations:
                    next(bg)
                    bgstr = BG.YELLOW
                elif (x, y) == move:
                    next(bg)
                    bgstr = BG.GREEN
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

    def print_board(self):
        os.system('cls' if os.name == 'nt' else 'clear')
        print(str(self))
        for line in self._log[-12:]:
            print(line)

    def __str__(self):
        ret = '   '
        for x in 'abcdefgh':
            ret += ' ' + x + ' '
        ret += '\n'
        bg = self.bg_iter()
        for y, row in enumerate(self.grid):
            ret += str(self.dim-y) + '  '
            for cell in row:
                if cell:
                    ret += cell.color.fg
                    ret += next(bg) + ' ' + cell.symbol + ' '
                else:
                    ret += next(bg) + '   '
            ret += CLEAR_FMT + '  ' + str(self.dim-y) + '\n'
        ret += '   '
        for x in 'abcdefgh':
            ret += ' ' + x + ' '
        ret += CLEAR_FMT
        return ret


class Piece:

    symbol = 'X'
    name = 'Generic Piece'

    def __init__(self, color):
        self.color = color
        self.dead = False

    def __str__(self):
        return self.name + ' (' + self.color.name + ')'
    
    def validate_move(self, x2, y2, board):
        x1, y1 = self.get_coords(board)
        try:
            assert((x1, y1) != (x2, y2))
            if board[x2,y2]:
                assert(board[x2,y2].color != self.color)
            return True
        except AssertionError:
            return False

    def get_coords(self, board):
        for i, piece in enumerate(board):
            if piece is self:
                return (i%board.dim, i//board.dim)
        
    def theoretical_moves(self):
        return []

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


class Rook(Piece):

    symbol = 'R'
    name = 'Rook'

    def validate_move(self, x2, y2, board):
        try:
            assert(Piece.validate_move(self, x2, y2, board))
            x1, y1 = self.get_coords(board)
            if y1 == y2:
                mvmt_op = op.add if x2 > x1 else op.sub
                x = mvmt_op(x1, 1)
                while x != x2:
                    assert(not board[x,y1])
                    x = mvmt_op(x, 1)
            elif x1 == x2:
                mvmt_op = op.add if y2 > y1 else op.sub
                y = mvmt_op(y1, 1)
                while y != y2:
                    assert(not board[x1,y])
                    y = mvmt_op(y, 1)
            else:
                return False
            return True
        except AssertionError:
            return False

    def theoretical_moves(self, board):
        x1, y1 = self.get_coords(board)
        for y in range(board.dim):
            for x in range(board.dim):
                if ((x == x1 or y == y1) and (x, y) != (x1, y1)):
                    yield (x, y)

class Bishop(Piece):

    symbol = 'B'
    name = 'Bishop'

    def validate_move(self, x2, y2, board):
        try:
            assert(Piece.validate_move(self, x2, y2, board))
            x1, y1 = self.get_coords(board)
            assert(abs(y2-y1) == abs(x2-x1))
            vert_op = op.add if y2 > y1 else op.sub 
            hor_op =  op.add if x2 > x1 else op.sub
            x,y = hor_op(x1, 1), vert_op(y1, 1)
            while (x,y) != (x2, y2):
                assert(not board[x,y])
                x,y = hor_op(x, 1), vert_op(y, 1)
            return True
        except AssertionError:
            return False
    def theoretical_moves(self, board):
        x1, y1 = self.get_coords(board)
        for y in range(board.dim):
            for x in range(board.dim):
                if ((abs(y-y1) == abs(x-x1)) and (x, y) != (x1, y1)):
                    yield (x, y)

class Queen(Piece):

    symbol = 'Q'
    name = 'Queen'

    def validate_move(self, x2, y2, board):
        try:
            x1, y1 = self.get_coords(board)
            assert(Piece.validate_move(self, x2, y2, board))
            if x1 == x2 or y1 == y2:
                assert(Rook.validate_move(self, x2, y2, board))
            elif abs(y2-y1) == abs(x2-x1):
                assert(Bishop.validate_move(self, x2, y2, board))
            else:
                return False
            return True
        except AssertionError:
            return False

    def theoretical_moves(self, board):
        x1, y1 = self.get_coords(board)
        for y in range(board.dim):
            for x in range(board.dim):
                if ((abs(y-y1) == abs(x-x1)) and (x, y) != (x1, y1)):
                    yield (x, y)
                if ((x == x1 or y == y1) and (x, y) != (x1, y1)):
                    yield (x, y)

class King(Piece):

    symbol = 'K'
    name = 'King'

    def validate_move(self, x2, y2, board):
        try:
            x1, y1 = self.get_coords(board)
            assert(Piece.validate_move(self, x2, y2, board))
            assert(abs(x2-x1) <= 1 and abs(y2-y1) <= 1)
            assert(Queen.validate_move(self, x2, y2, board))
            return True
        except AssertionError:
            return False

    def theoretical_moves(self, board):
        x, y = self.get_coords(board)
        yield (x+1, y)
        yield (x+1, y+1)
        yield (x, y+1)
        yield (x-1, y+1)
        yield (x-1, y)
        yield (x-1, y-1)
        yield (x, y-1)
        yield (x+1, y-1)

class Knight(Piece):

    symbol = 'N'
    name = 'Knight'

    def validate_move(self, x2, y2, board):
        try:
            x1, y1 = self.get_coords(board)
            assert(Piece.validate_move(self, x2, y2, board))
            assert((abs(y2-y1) == 2 and abs(x2-x1) == 1) or (abs(y2-y1) == 1 and abs(x2-x1) == 2))
            return True
        except AssertionError:
            return False

    def theoretical_moves(self, board):
        x1, y1 = self.get_coords(board)
        for y in range(board.dim):
            for x in range(board.dim):
                if ((abs(y-y1) == 2 and abs(x-x1) == 1) or (abs(y-y1) == 1 and abs(x-x1) == 2)):
                    yield (x, y)


inp = input('AI? y/[n]')
ai = inp.lower() == 'y'
board = ChessBoard(AI=ai)
board.play()

