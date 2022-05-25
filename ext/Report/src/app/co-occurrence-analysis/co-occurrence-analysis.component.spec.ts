import {ComponentFixture, TestBed} from '@angular/core/testing';

import {CoOccurrenceAnalysisComponent} from './co-occurrence-analysis.component';

describe('CoOccurrenceAnalysusComponent', () => {
  let component: CoOccurrenceAnalysisComponent;
  let fixture: ComponentFixture<CoOccurrenceAnalysisComponent>;

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      declarations: [CoOccurrenceAnalysisComponent]
    })
      .compileComponents();
  });

  beforeEach(() => {
    fixture = TestBed.createComponent(CoOccurrenceAnalysisComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
