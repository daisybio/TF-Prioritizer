import {NgModule} from '@angular/core';
import {BrowserModule} from '@angular/platform-browser';

import {AppComponent} from './app.component';
import {BrowserAnimationsModule} from '@angular/platform-browser/animations';
import {MatButtonModule} from "@angular/material/button";
import {RouterOutlet} from "@angular/router";
import {AppRoutingModule} from './app-routing.module';
import {RankingComponent} from "./pages/ranking/ranking.component";
import {MatAutocompleteModule} from "@angular/material/autocomplete";
import {MatFormFieldModule} from "@angular/material/form-field";
import {MatInputModule} from "@angular/material/input";
import {MatChipsModule} from "@angular/material/chips";
import {MatExpansionModule} from "@angular/material/expansion";
import {MatDividerModule} from "@angular/material/divider";
import {MatTableModule} from "@angular/material/table";
import {BarChartModule, HeatMapModule} from "@swimlane/ngx-charts";
import {FormsModule, ReactiveFormsModule} from "@angular/forms";
import {HashLocationStrategy, LocationStrategy} from "@angular/common";
import {CoOccurrenceComponent} from './pages/generalAnalysis/co-occurrence/co-occurrence.component';
import {
  StatisticalEvaluationComponent
} from './pages/generalAnalysis/statistical-evaluation/statistical-evaluation.component';
import {ParametersComponent} from './pages/explainatory/parameters/parameters.component';
import {DocumentationComponent} from './pages/explainatory/documentation/documentation.component';
import {ValidationComponent} from './pages/perTf/validation/validation.component';
import {DistributionComponent} from './pages/perTf/distribution/distribution.component';
import {RegressionComponent} from './pages/perTf/regression/regression.component';
import {MatTooltipModule} from "@angular/material/tooltip";
import {ImageSelectorComponent} from './components/image-selector/image-selector.component';
import {MatMenuModule} from "@angular/material/menu";
import {MatSelectModule} from "@angular/material/select";
import {TfDetailsComponent} from './pages/perTf/tf-details/tf-details.component';
import {MatTabsModule} from "@angular/material/tabs";
import {GeneralInformationComponent} from './components/general-information/general-information.component';
import {
  RegressionCoefficientsComponent
} from './pages/generalAnalysis/regression-coefficients/regression-coefficients.component';
import {DataSelectorComponent} from './components/data-selector/data-selector.component';
import {MatSliderModule} from "@angular/material/slider";
import {MatListModule} from "@angular/material/list";
import {MatButtonToggleModule} from "@angular/material/button-toggle";
import {MatIconModule} from "@angular/material/icon";

@NgModule({
  declarations: [
    AppComponent,
    RankingComponent,
    CoOccurrenceComponent,
    StatisticalEvaluationComponent,
    ParametersComponent,
    DocumentationComponent,
    ValidationComponent,
    DistributionComponent,
    RegressionComponent,
    ImageSelectorComponent,
    TfDetailsComponent,
    GeneralInformationComponent,
    RegressionCoefficientsComponent,
    DataSelectorComponent,
  ],
  imports: [
    BrowserModule,
    BrowserAnimationsModule,
    MatButtonModule,
    RouterOutlet,
    AppRoutingModule,
    MatAutocompleteModule,
    MatFormFieldModule,
    MatInputModule,
    MatChipsModule,
    MatExpansionModule,
    MatDividerModule,
    MatTableModule,
    BarChartModule,
    ReactiveFormsModule,
    MatTooltipModule,
    MatMenuModule,
    MatSelectModule,
    MatTabsModule,
    MatSliderModule,
    FormsModule,
    MatListModule,
    HeatMapModule,
    MatButtonToggleModule,
    MatIconModule,
  ],
  providers: [
    {
      provide: LocationStrategy, useClass: HashLocationStrategy
    }
  ],
  bootstrap: [AppComponent]
})
export class AppModule {
}
