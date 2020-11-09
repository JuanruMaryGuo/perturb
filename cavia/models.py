import torch
import torch.nn as nn
import torch.nn.functional as F


class MLP(nn.Module):

    def __init__(self,
                 num_classes = 10,
                 input_size=100,
                 num_context_params = 100,
                 num_context_params_layers = 0,
                 num_film_hidden_layers = 0,
                 initialisation = 'xavier',
                 context_in = [False, False, False, False, False],
                 hidden_size1=80, hidden_size2=100, device = "cpu"):

        super(MLP, self).__init__()
        self.num_classes = num_classes
        self.input_size = input_size
        self.num_context_params = num_context_params
        self.num_film_hidden_layers = num_film_hidden_layers
        self.num_context_params_laysers = num_context_params_layers
        self.hidden_size1 = hidden_size1
        self.hidden_size2 = hidden_size2
        self.context_in = context_in
        self.mlp1 = nn.Linear(self.input_size, self.hidden_size1).to(device)
        self.mlp2 = nn.Linear(self.hidden_size1, self.hidden_size2).to(device)
        self.mlp3 = nn.Linear(self.hidden_size2, self.hidden_size1).to(device)
        self.mlp4 = nn.Linear(self.hidden_size1, self.hidden_size2).to(device)
        self.mlp5 = nn.Linear(self.hidden_size2, self.hidden_size1).to(device)
        self.classify = nn.Linear(self.hidden_size1,self.num_classes).to(device)

        # for each layer where we have context parameters, initialise a FiLM layer
        if self.context_in[0]:
            self.film1 = nn.Linear(self.num_context_params, 2 * self.hidden_size1).to(device)
            if self.num_film_hidden_layers == 1:
                self.film11 = nn.Linear(2, 2).to(device)
        if self.context_in[1]:
            self.film2 = nn.Linear(self.num_context_params, 2 * self.hidden_size2).to(device)
            if self.num_film_hidden_layers == 1:
                self.film22 = nn.Linear(2, 2).to(device)
        if self.context_in[2]:
            self.film3 = nn.Linear(self.num_context_params, 2 * self.hidden_size1).to(device)
            if self.num_film_hidden_layers == 1:
                self.film33 = nn.Linear(2, 2).to(device)
        if self.context_in[3]:
            self.film4 = nn.Linear(self.num_context_params, 2 * self.hidden_size2).to(device)
            if self.num_film_hidden_layers == 1:
                self.film44 = nn.Linear(2, 2).to(device)


        # parameter initialisation (if different than standard pytorch one)
        if initialisation != 'standard':
            self.init_params(initialisation)

        # initialise context parameters
        self.context_params = torch.zeros(size=[self.num_context_params], requires_grad=True).to(device)

    def init_params(self, initialisation):


        if initialisation == 'xavier':
            torch.nn.init.xavier_uniform_(self.mlp1.weight, gain=nn.init.calculate_gain('relu', self.mlp1.weight))
            torch.nn.init.xavier_uniform_(self.mlp2.weight, gain=nn.init.calculate_gain('relu', self.mlp2.weight))
            torch.nn.init.xavier_uniform_(self.mlp3.weight, gain=nn.init.calculate_gain('relu', self.mlp3.weight))
            torch.nn.init.xavier_uniform_(self.mlp4.weight, gain=nn.init.calculate_gain('relu', self.mlp4.weight))
            torch.nn.init.xavier_uniform_(self.mlp5.weight, gain=nn.init.calculate_gain('relu', self.mlp5.weight))
        elif initialisation == 'kaiming':
            torch.nn.init.kaiming_uniform_(self.mlp1.weight, nonlinearity='relu')
            torch.nn.init.kaiming_uniform_(self.mlp2.weight, nonlinearity='relu')
            torch.nn.init.kaiming_uniform_(self.mlp3.weight, nonlinearity='relu')
            torch.nn.init.kaiming_uniform_(self.mlp4.weight, nonlinearity='relu')
            torch.nn.init.kaiming_uniform_(self.mlp5.weight, nonlinearity='relu')


        self.mlp1.bias.data.fill_(0)
        self.mlp2.bias.data.fill_(0)
        self.mlp3.bias.data.fill_(0)
        self.mlp4.bias.data.fill_(0)
        self.mlp5.bias.data.fill_(0)

        # fully connected weights at the end

        if initialisation == 'xavier':
            torch.nn.init.xavier_uniform_(self.classify.weight, gain=nn.init.calculate_gain('linear', self.classify.weight))
        elif initialisation == 'kaiming':
            torch.nn.init.kaiming_uniform_(self.classify.weight, nonlinearity='linear')

        # fully connected bias

        self.classify.bias.data.fill_(0)

        # FiLM layer weights

        if self.context_in[0] and initialisation == 'xavier':
            torch.nn.init.xavier_uniform_(self.film1.weight, gain=nn.init.calculate_gain('linear', self.film1.weight))
        elif self.context_in[0] and initialisation == 'kaiming':
            torch.nn.init.kaiming_uniform_(self.film1.weight, nonlinearity='linear')

        if self.context_in[1] and initialisation == 'xavier':
            torch.nn.init.xavier_uniform_(self.film2.weight, gain=nn.init.calculate_gain('linear', self.film2.weight))
        elif self.context_in[1] and initialisation == 'kaiming':
            torch.nn.init.kaiming_uniform_(self.film2.weight, nonlinearity='linear')

        if self.context_in[2] and initialisation == 'xavier':
            torch.nn.init.xavier_uniform_(self.film3.weight, gain=nn.init.calculate_gain('linear', self.film3.weight))
        elif self.context_in[2] and initialisation == 'kaiming':
            torch.nn.init.kaiming_uniform_(self.film3.weight, nonlinearity='linear')

        if self.context_in[3] and initialisation == 'xavier':
            torch.nn.init.xavier_uniform_(self.film4.weight, gain=nn.init.calculate_gain('linear', self.film4.weight))
        elif self.context_in[3] and initialisation == 'kaiming':
            torch.nn.init.kaiming_uniform_(self.film4.weight, nonlinearity='linear')


    def reset_context_params(self):
        self.context_params = self.context_params.detach() * 0
        self.context_params.requires_grad = True

    def forward(self, x):

        h1 = self.mlp1(x)
        if self.context_in[0]:
            film1 = self.film1(self.context_params)
            if self.num_context_params_laysers == 1:
                film1 = self.film11(F.relu(film1))
            gamma1 = film1[:self.hidden_size1]
            beta1 = film1[self.hidden_size1:]
            h1 = gamma1 * h1 + beta1
        h1 = F.relu(h1)

        h2 = self.mlp2(h1)
        if self.context_in[1]:
            film2 = self.film2(self.context_params)
            if self.num_film_hidden_layers == 1:
                film2 = self.film22(F.relu(film2))
            gamma2 = film2[:self.hidden_size2]
            beta2 = film2[self.hidden_size2:]
            h2 = gamma2 * h2 + beta2
        h2 = F.relu(h2)

        h3 = self.mlp3(h2)
        if self.context_in[2]:
            film3 = self.film3(self.context_params)
            if self.num_film_hidden_layers == 1:
                film3 = self.film33(F.relu(film3))
            gamma3 = film3[:self.hidden_size1]
            beta3 = film3[self.hidden_size1:]
            h3 = gamma3 * h3 + beta3
        h3 = F.relu(h3)

        h4 = self.mlp4(h3)
        if self.context_in[3]:
            film4 = self.film4(self.context_params)
            if self.num_film_hidden_layers == 1:
                film4 = self.film44(F.relu(film4))
            gamma4 = film4[:self.hidden_size2]
            beta4 = film4[self.hidden_size2:]
            h4 = gamma4 * h4 + beta4
        h4 = F.relu(h4)

        h5 = self.mlp5(h4)
        h5 = F.relu(h5)

        y = self.classify(h5)

        return y

# for test
if __name__ == '__main__':
    net = CondConvNet()
    print(net)
    print(list(net.parameters()))
